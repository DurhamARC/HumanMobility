#!/usr/bin/env python3
"""
Performance Analysis Visualization (pa_vis.py)

Creates a two-page PDF:
- Page 1: Memory Usage and CPU Usage across MPI ranks
- Page 2: Thread Usage with bars for every individual thread

Usage:
    python pa_vis.py [--use-balance-logs] [--use-dat-files] [--use-thread-files]
    
By default, tries to use balance logs first, then falls back to .dat files.
"""

import numpy as np
import sys
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def parse_balance_logs():
    """Parse rank_memory_balance.log and rank_cpu_balance.log files"""
    memfile = "rank_memory_balance.log"
    cpufile = "rank_cpu_balance.log"
    
    if not (os.path.exists(memfile) and os.path.exists(cpufile)):
        return None, None
    
    print(f"Using balance log files: {memfile}, {cpufile}")
    
    # Load memory data
    memdat = np.loadtxt(
        memfile,
        dtype=[("step", int), ("rank", int), ("resident", float)],
    )
    
    # Load CPU data  
    cpudat = np.loadtxt(
        cpufile,
        dtype=[("step", int), ("rank", int), ("user", float), 
               ("sys", float), ("sum", float), ("deadfrac", float)],
    )
    
    return memdat, cpudat

def parse_thread_info_files():
    """Parse thread_info_MPI-step*.dat files for thread usage analysis"""
    thread_files = glob.glob("thread_info_MPI-step*.dat")
    
    if not thread_files:
        print("No thread_info_MPI files found")
        return None
    
    print(f"Using thread info files: {len(thread_files)} files found")
    
    # Parse thread timing data
    thread_data = {}  # rank -> {step -> {thread_id -> total_time}}
    rank_thread_counts = {}  # rank -> max_threads
    
    for thread_file in sorted(thread_files):
        # Extract step from filename: thread_info_MPI-step1.dat
        step_part = thread_file.split('-step')[1].split('.dat')[0]
        step = int(step_part)
        
        try:
            # Load thread data - format from SWIFT source code analysis
            data = np.loadtxt(thread_file)
            
            if data.size == 0:
                continue
                
            # First row contains metadata: rank, tic_step, toc_step, etc.
            if len(data.shape) == 1:
                # Single row file
                continue
                
            # Skip first row (metadata) and process task data
            task_data = data[1:] if data.shape[0] > 1 else data
            
            if task_data.size == 0:
                continue
            
            # Process each task line
            # Format: rank rid type subtype pair tic toc ci.hydro.count cj.hydro.count ci.grav.count cj.grav.count flags sid
            if len(task_data.shape) == 2:
                for row in task_data:
                    if len(row) >= 7:  # Ensure we have enough columns
                        rank = int(row[0])  # First column is rank
                        thread_id = int(row[1])  # rid = runner/thread ID
                        tic = int(row[5])
                        toc = int(row[6])
                        
                        if toc > tic:  # Valid task timing
                            task_time = toc - tic  # Time in ticks
                            
                            # Initialize rank data structure
                            if rank not in thread_data:
                                thread_data[rank] = {}
                                rank_thread_counts[rank] = 0
                            
                            if step not in thread_data[rank]:
                                thread_data[rank][step] = {}
                            
                            if thread_id not in thread_data[rank][step]:
                                thread_data[rank][step][thread_id] = 0
                            
                            thread_data[rank][step][thread_id] += task_time
                            rank_thread_counts[rank] = max(rank_thread_counts[rank], thread_id + 1)
            
        except Exception as e:
            print(f"Warning: Could not parse {thread_file}: {e}")
            continue
    
    if not thread_data:
        return None
    
    return thread_data, rank_thread_counts

def parse_dat_files():
    """Parse memuse_report-rank*-step*.dat and mpiuse_report-rank*-step*.dat files"""
    memfiles = glob.glob("memuse_report-rank*-step*.dat")
    cpufiles = glob.glob("mpiuse_report-rank*-step*.dat")  # MPI communication files
    
    if not memfiles:
        print("No memuse_report files found")
        return None, None
    
    print(f"Using .dat files: {len(memfiles)} memuse files, {len(cpufiles)} mpiuse files")
    
    # Parse memory usage files
    mem_data = []
    for memfile in sorted(memfiles):
        # Extract rank and step from filename: memuse_report-rank0-step1.dat
        parts = os.path.basename(memfile).replace('.dat', '').split('-')
        rank = int(parts[1].replace('rank', ''))
        step = int(parts[2].replace('step', ''))
        
        try:
            # Load memory data - format varies, so we need to be careful
            with open(memfile, 'r') as f:
                lines = f.readlines()
                
            # Find peak memory usage from comments or data
            peak_mem = 0
            for line in lines:
                if line.startswith('# Peak memory usage'):
                    # Extract MB value and convert to KB
                    peak_mb = float(line.split(':')[1].strip().split()[0])
                    peak_mem = peak_mb * 1024  # Convert MB to KB
                    break
            
            if peak_mem == 0:
                # Try to get from process memory line
                for line in lines:
                    if 'Memory use by process' in line:
                        # Try to extract memory value - this is system dependent
                        try:
                            parts = line.split()
                            for i, part in enumerate(parts):
                                if 'MB' in part or 'KB' in part:
                                    val = float(parts[i-1])
                                    if 'MB' in part:
                                        peak_mem = val * 1024
                                    else:
                                        peak_mem = val
                                    break
                        except:
                            peak_mem = 1000000  # Default 1GB in KB
                        break
            
            if peak_mem == 0:
                peak_mem = 1000000  # Default 1GB in KB
                
            mem_data.append((step, rank, peak_mem))
            
        except Exception as e:
            print(f"Warning: Could not parse {memfile}: {e}")
            continue
    
    # Parse MPI/CPU usage files (these represent communication time, not total CPU)
    cpu_data = []
    for cpufile in sorted(cpufiles):
        # Extract rank and step from filename
        parts = os.path.basename(cpufile).replace('.dat', '').split('-')
        rank = int(parts[1].replace('rank', ''))
        step = int(parts[2].replace('step', ''))
        
        try:
            # Load MPI communication data
            with open(cpufile, 'r') as f:
                lines = f.readlines()
            
            # Sum up communication times as a proxy for CPU usage
            total_comm_time = 0
            for line in lines:
                if line.startswith('#') or not line.strip():
                    continue
                try:
                    # Format: stic etic dtic step rank otherrank type itype subtype isubtype activation tag size sum
                    parts = line.split()
                    if len(parts) >= 3:
                        dtic = float(parts[2])  # Duration of MPI operation
                        total_comm_time += dtic
                except:
                    continue
            
            # Convert ticks to milliseconds (approximate)
            cpu_time_ms = total_comm_time / 1000.0  # Rough conversion
            cpu_data.append((step, rank, 0, 0, cpu_time_ms, 0))  # user, sys, sum, deadfrac
            
        except Exception as e:
            print(f"Warning: Could not parse {cpufile}: {e}")
            continue
    
    if not mem_data and not cpu_data:
        return None, None
    
    # Convert to numpy arrays with same format as balance logs
    if mem_data:
        memdat = np.array(mem_data, dtype=[("step", int), ("rank", int), ("resident", float)])
    else:
        memdat = None
        
    if cpu_data:
        cpudat = np.array(cpu_data, dtype=[
            ("step", int), ("rank", int), ("user", float), 
            ("sys", float), ("sum", float), ("deadfrac", float)
        ])
    else:
        cpudat = None
    
    return memdat, cpudat

def get_node_assignments(num_ranks):
    """Automatically assign node categories based on number of ranks"""
    if num_ranks <= 8:
        return ["background"] * num_ranks
    else:
        mid_point = num_ranks // 2
        return ["background"] * mid_point + ["zoom"] * (num_ranks - mid_point)

def main():
    parser = argparse.ArgumentParser(description='Performance Analysis Visualization')
    parser.add_argument('--use-balance-logs', action='store_true', 
                       help='Force use of rank_*_balance.log files')
    parser.add_argument('--use-dat-files', action='store_true',
                       help='Force use of *_report-rank*-step*.dat files')
    parser.add_argument('--use-thread-files', action='store_true',
                       help='Force use of thread_info_MPI-step*.dat files')
    args = parser.parse_args()
    
    # Color mapping
    cmap = {"zoom": "tab:red", "background": "tab:blue", "compute": "tab:green"}
    
    # Try to load data
    memdat, cpudat = None, None
    thread_data, rank_thread_counts = None, None
    
    if args.use_thread_files:
        thread_result = parse_thread_info_files()
        if thread_result:
            thread_data, rank_thread_counts = thread_result
    elif args.use_dat_files:
        memdat, cpudat = parse_dat_files()
        thread_result = parse_thread_info_files()
        if thread_result:
            thread_data, rank_thread_counts = thread_result
    elif args.use_balance_logs:
        memdat, cpudat = parse_balance_logs()
    else:
        # Try balance logs first, then .dat files, then thread files
        memdat, cpudat = parse_balance_logs()
        if memdat is None or cpudat is None:
            print("Balance logs not found, trying .dat files...")
            memdat, cpudat = parse_dat_files()
        
        # Always try to get thread data
        thread_result = parse_thread_info_files()
        if thread_result:
            thread_data, rank_thread_counts = thread_result
    
    if memdat is None and cpudat is None and thread_data is None:
        print("Error: No suitable data files found!")
        print("Looking for:")
        print("  - rank_memory_balance.log and rank_cpu_balance.log")
        print("  - memuse_report-rank*-step*.dat files")
        print("  - mpiuse_report-rank*-step*.dat files")
        print("  - thread_info_MPI-step*.dat files")
        sys.exit(1)
    
    # Get unique ranks and steps
    if memdat is not None:
        unique_ranks = np.unique(memdat["rank"])
        memory_steps = np.unique(memdat["step"])
    else:
        memory_steps = np.array([])
        
    if cpudat is not None:
        if memdat is None:
            unique_ranks = np.unique(cpudat["rank"])
        cpu_steps = np.unique(cpudat["step"])
    else:
        cpu_steps = np.array([])
        
    if thread_data is not None:
        if memdat is None and cpudat is None:
            unique_ranks = np.array(sorted(thread_data.keys()))
        all_thread_steps = set()
        for rank_steps in thread_data.values():
            all_thread_steps.update(rank_steps.keys())
        thread_steps = np.array(sorted(all_thread_steps))
    else:
        thread_steps = np.array([])
    
    # Use appropriate step count for each page
    balance_steps = np.union1d(memory_steps, cpu_steps)
    
    num_ranks = len(unique_ranks)
    num_balance_steps = len(balance_steps)
    num_thread_steps = len(thread_steps)
    
    print(f"Found {num_ranks} ranks and {num_balance_steps} balance steps")
    print(f"Ranks: {unique_ranks}")
    print(f"Balance Steps: {balance_steps}")
    
    # Set up node assignments
    node_assignments = get_node_assignments(num_ranks)
    
    # Calculate accumulated values per rank
    accumulated_memory = np.zeros(num_ranks)
    accumulated_cpu = np.zeros(num_ranks)
    accumulated_thread_time = np.zeros(num_ranks)
    thread_efficiency = np.zeros(num_ranks)
    max_memory = np.zeros(num_ranks)
    
    # Individual thread data for detailed visualization
    individual_thread_data = {}  # rank -> {thread_id -> total_time}
    
    for i, rank in enumerate(unique_ranks):
        # Memory analysis
        if memdat is not None:
            rank_mask_mem = memdat["rank"] == rank
            rank_memory = memdat["resident"][rank_mask_mem]
            if len(rank_memory) > 0:
                accumulated_memory[i] = np.mean(rank_memory)  # Average memory usage
                max_memory[i] = np.max(rank_memory)  # Peak memory usage
        
        # CPU analysis
        if cpudat is not None:
            rank_mask_cpu = cpudat["rank"] == rank
            rank_cpu = cpudat["sum"][rank_mask_cpu]
            if len(rank_cpu) > 0:
                accumulated_cpu[i] = np.sum(rank_cpu)  # Total CPU time
        
        # Thread analysis
        if thread_data is not None and rank in thread_data:
            total_thread_time = 0
            max_threads = rank_thread_counts.get(rank, 1)
            thread_utilization = 0
            
            # Initialize individual thread data for this rank
            individual_thread_data[rank] = {}
            
            for step, thread_times in thread_data[rank].items():
                step_total = sum(thread_times.values())
                total_thread_time += step_total
                
                # Accumulate individual thread times
                for thread_id, thread_time in thread_times.items():
                    if thread_id not in individual_thread_data[rank]:
                        individual_thread_data[rank][thread_id] = 0
                    individual_thread_data[rank][thread_id] += thread_time
                
                # Calculate thread utilization for this step
                if thread_times:
                    active_threads = len(thread_times)
                    thread_utilization += active_threads / max_threads
            
            accumulated_thread_time[i] = total_thread_time / 1000000000.0  # Convert ticks to seconds
            if len(thread_data[rank]) > 0:
                thread_efficiency[i] = thread_utilization / len(thread_data[rank])  # Average thread efficiency
    
    # Set up plotting labels and colors
    bar_labels = []
    bar_colors = []
    for i, node_assignment in enumerate(node_assignments):
        if node_assignment in bar_labels:
            bar_labels.append(f"_{node_assignment}")
        else:
            bar_labels.append(node_assignment)
        bar_colors.append(cmap[node_assignment])
    
    # Create the two-page visualization
    print("Creating performance analysis visualization...")
    
    with PdfPages("pa_vis.pdf") as pdffile:
        
        # PAGE 1: Memory Usage and CPU Usage
        has_memory_data = memdat is not None and np.any(accumulated_memory > 0)
        has_cpu_data = cpudat is not None and np.any(accumulated_cpu > 0)

        if has_memory_data or has_cpu_data:
            # Calculate number of plots needed
            num_plots = int(has_memory_data) + int(has_cpu_data)
            
            # Fix: Increase figure width to accommodate both charts properly
            fig1, axs1 = plt.subplots(1, 2, figsize=(16, 6))
            
            # Fix: Center the title properly by adjusting spacing and using figure-level suptitle
            fig1.suptitle(f"Memory and CPU Performance Analysis\n{num_balance_steps} selected steps across {num_ranks} MPI ranks", 
                         fontsize=14, ha='center', va='top', y=0.95)
            
            plot_idx = 0
            
            # Memory usage chart
            if has_memory_data:
                bars_mem = axs1[plot_idx].bar(unique_ranks, accumulated_memory, color=bar_colors, alpha=0.7)
                axs1[plot_idx].axhline(np.mean(accumulated_memory), color="black", ls="dashed", lw=1.5, label="Average")
                axs1[plot_idx].axhline(1024**2, color="red", ls="solid", lw=2, label="1 GB threshold") # 1024^2 KB
                axs1[plot_idx].set_ylabel("Average Memory Usage [KB]")
                axs1[plot_idx].set_xlabel("MPI Rank")
                axs1[plot_idx].set_title("Memory Usage (Average over selected steps)")
                axs1[plot_idx].legend(loc="lower right", fontsize=9)
                axs1[plot_idx].grid(True, alpha=0.3)
                
                # Add text annotations for memory values
                for i, (rank, mem_kb) in enumerate(zip(unique_ranks, accumulated_memory)):
                    mem_gb = mem_kb / (1024**2)
                    axs1[plot_idx].text(rank, mem_kb + 0.08*(1024**2), f"{mem_gb:.2f}GB", 
                               ha='center', va='bottom', fontsize=10)

                # Set Y-axis limits to provide more space
                max_mem = np.max(accumulated_memory)
                axs1[plot_idx].set_ylim(0, max_mem * 1.15)  # Add 15% padding at top
                axs1[plot_idx].set_xticks(unique_ranks)
                plot_idx += 1
            
            # CPU usage chart
            if has_cpu_data:
                bars_cpu = axs1[plot_idx].bar(unique_ranks, accumulated_cpu / 1000,  # Convert to seconds
                                     color=bar_colors, alpha=0.7)
                axs1[plot_idx].axhline(np.mean(accumulated_cpu) / 1000, color="black", 
                               ls="dashed", lw=1.5, label="Average")
                axs1[plot_idx].set_ylabel("Total CPU Time [seconds]")
                axs1[plot_idx].set_xlabel("MPI Rank")
                axs1[plot_idx].set_title("CPU Usage (Total over selected steps)")
                
                # Calculate and display efficiency
                total_cpu_time = accumulated_cpu.sum()
                max_cpu_time = accumulated_cpu.max()
                efficiency = total_cpu_time / (num_ranks * max_cpu_time) if max_cpu_time > 0 else 0
                
                axs1[plot_idx].text(0.95, 0.95, 
                           f"ε = {100*efficiency:.1f}%\nTotal: {total_cpu_time/1000:.1f}s",
                           transform=axs1[plot_idx].transAxes, va="top", ha="right",
                           bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.7))
                
                axs1[plot_idx].legend(loc="lower right", fontsize=9)
                axs1[plot_idx].grid(True, alpha=0.3)
                
                # Add text annotations for CPU values
                for i, (rank, cpu_s) in enumerate(zip(unique_ranks, accumulated_cpu / 1000)):
                    axs1[plot_idx].text(rank, cpu_s + max(accumulated_cpu)/1000 * 0.05, f"{cpu_s:.1f}", 
                               ha='center', va='bottom', fontsize=10)

                # Set Y-axis limits to provide more space
                max_cpu_s = np.max(accumulated_cpu) / 1000
                axs1[plot_idx].set_ylim(0, max_cpu_s * 1.15)  # Add 15% padding at top
                axs1[plot_idx].set_xticks(unique_ranks)
                plot_idx += 1
        
        # Hide unused subplot if only one plot
        if plot_idx == 1:
            axs1[1].set_visible(False)
        
        # Fix: Adjust spacing for wider figure and better layout
        plt.subplots_adjust(top=0.85, bottom=0.15, left=0.08, right=0.95, wspace=0.25)
        plt.savefig(pdffile, format="pdf", bbox_inches='tight', dpi=150)
        plt.close()
        
        # PAGE 2: Individual Thread Usage (separated by rank)
        if thread_data is not None and individual_thread_data:
            # Determine subplot layout based on number of ranks
            num_ranks_with_threads = len(individual_thread_data)
            
            if num_ranks_with_threads == 1:
                fig2, ax2 = plt.subplots(1, 1, figsize=(12, 6))
                axs2 = [ax2]
            elif num_ranks_with_threads == 2:
                fig2, axs2 = plt.subplots(1, 2, figsize=(16, 6))
            elif num_ranks_with_threads <= 4:
                fig2, axs2 = plt.subplots(2, 2, figsize=(16, 10))
                axs2 = axs2.flatten()
            else:
                # For more than 4 ranks, create a grid layout
                ncols = min(3, num_ranks_with_threads)
                nrows = (num_ranks_with_threads + ncols - 1) // ncols
                fig2, axs2 = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
                if nrows > 1:
                    axs2 = axs2.flatten()
                else:
                    axs2 = [axs2] if ncols == 1 else axs2
            
            fig2.suptitle(f"Individual Thread Usage by Rank\n{num_thread_steps} selected steps across {num_ranks} MPI ranks", 
                         fontsize=14)
            
            # Plot thread usage for each rank separately
            rank_idx = 0
            global_max_thread_time = 0
            
            # First pass: find global maximum for consistent y-axis scaling
            for rank in sorted(individual_thread_data.keys()):
                for thread_time in individual_thread_data[rank].values():
                    thread_time_s = thread_time / 1000000000.0  # Convert to seconds
                    global_max_thread_time = max(global_max_thread_time, thread_time_s)
            
            # Second pass: create plots
            for rank in sorted(individual_thread_data.keys()):
                if rank_idx >= len(axs2):
                    break
                    
                ax = axs2[rank_idx]
                
                # Get thread data for this rank
                rank_threads = individual_thread_data[rank]
                thread_ids = sorted(rank_threads.keys())
                thread_times = [rank_threads[tid] / 1000000000.0 for tid in thread_ids]  # Convert to seconds
                thread_labels = [f"T{tid}" for tid in thread_ids]
                
                if thread_times:
                    # Create bar chart for this rank's threads
                    x_positions = range(len(thread_ids))
                    rank_color = bar_colors[rank] if rank < len(bar_colors) else 'tab:gray'
                    bars = ax.bar(x_positions, thread_times, color=rank_color, alpha=0.7)
                    
                    # Add average line for this rank
                    rank_avg = np.mean(thread_times)
                    ax.axhline(rank_avg, color="black", ls="dashed", lw=1.5, label=f"Rank {rank} Avg")
                    
                    # Set labels and title
                    ax.set_ylabel("Thread Time [s]")
                    ax.set_xlabel("Thread ID")
                    ax.set_title(f"Rank {rank} Thread Usage")
                    ax.legend(loc="lower right", fontsize=9)
                    ax.grid(True, alpha=0.3)
                    
                    # Set x-axis labels
                    ax.set_xticks(x_positions)
                    ax.set_xticklabels(thread_labels)
                    
                    # Add text annotations for thread values
                    for i, (x_pos, thread_time) in enumerate(zip(x_positions, thread_times)):
                        if thread_time > 0:  # Only label non-zero values
                            ax.text(x_pos, thread_time + global_max_thread_time * 0.02, f"{thread_time:.1f}", 
                                   ha='center', va='bottom', fontsize=8, rotation=90)

                    # Calculate thread balance metrics for this rank
                    rank_total_time = sum(thread_times)
                    rank_max_time = max(thread_times) if thread_times else 0
                    rank_min_time = min(thread_times) if thread_times else 0
                    rank_balance_ratio = rank_min_time / rank_max_time if rank_max_time > 0 else 0
                    
                    # Add efficiency information
                    ax.text(0.98, 0.98, 
                            f"Total: {rank_total_time:.1f} s\n"
                            f"Balance: {rank_balance_ratio:.2f}\n"
                            f"Threads: {len(thread_times)}",
                            transform=ax.transAxes, va="top", ha="right",
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.7),
                            fontsize=9)
                    
                    # Set consistent Y-axis limits across all subplots
                    ax.set_ylim(0, global_max_thread_time * 1.15)
                else:
                    # No thread data for this rank
                    ax.text(0.5, 0.5, f"No thread data\nfor Rank {rank}", 
                           transform=ax.transAxes, ha='center', va='center',
                           bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.7))
                    ax.set_title(f"Rank {rank} Thread Usage")
                    ax.set_xlabel("Thread ID")
                    ax.set_ylabel("Thread Time [s]")
                
                rank_idx += 1
            
            # Hide unused subplots
            for i in range(rank_idx, len(axs2)):
                axs2[i].set_visible(False)
            
            plt.tight_layout()
            plt.savefig(pdffile, format="pdf", bbox_inches='tight', dpi=150)
            plt.close()
    
    # Print summary statistics
    print("\n" + "="*60)
    print("PERFORMANCE ANALYSIS SUMMARY")
    print("="*60)
    print(f"Simulation: {num_balance_steps} steps, {num_ranks} MPI ranks")

    # Show full list for balance steps, first/last for thread steps
    if num_balance_steps > 0:
        balance_steps_list = ", ".join(map(str, balance_steps))
        print(f"Memory & CPU steps analyzed: {balance_steps_list}")
    else:
        print("Memory & CPU steps analyzed: None")

    if num_thread_steps > 0:
        if num_thread_steps <= 10:
            # Show all steps if not too many
            thread_steps_list = ", ".join(map(str, thread_steps))
            print(f"Thread steps analyzed: {thread_steps_list}")
        else:
            # Show first and last few steps
            print(f"Thread steps analyzed: {thread_steps[0]} to {thread_steps[-1]} ({num_thread_steps} total steps)")
            print(f"  First 5: {', '.join(map(str, thread_steps[:5]))}")
            print(f"  Last 5: {', '.join(map(str, thread_steps[-5:]))}")
    else:
        print("Thread steps analyzed: None")

    print()
    
    if memdat is not None and np.any(accumulated_memory > 0):
        print("MEMORY USAGE:")
        for i, rank in enumerate(unique_ranks):
            print(f"  Rank {rank}: Avg = {accumulated_memory[i]/1024:.2f} GB, "
                  f"Max = {max_memory[i]/1024:.2f} GB")
        print(f"  Overall average: {np.mean(accumulated_memory)/1024:.2f} GB")
        print(f"  Memory balance: {np.std(accumulated_memory)/np.mean(accumulated_memory):.1%}")
    else:
        print("MEMORY USAGE: No data available")
    
    print()
    
    if cpudat is not None and np.any(accumulated_cpu > 0):
        print("CPU USAGE:")
        total_cpu_time = accumulated_cpu.sum()
        max_cpu_time = accumulated_cpu.max()
        efficiency = total_cpu_time / (num_ranks * max_cpu_time) if max_cpu_time > 0 else 0
        
        for i, rank in enumerate(unique_ranks):
            print(f"  Rank {rank}: Total = {accumulated_cpu[i]/1000:.1f} seconds")
        print(f"  Overall total: {total_cpu_time/1000:.1f} seconds")
        print(f"  Parallel efficiency: {100*efficiency:.1f}%")
        if max_cpu_time > 0:
            print(f"  Load balance ratio: {accumulated_cpu.min()/max_cpu_time:.2f}")
    else:
        print("CPU USAGE: No data available")
    
    print()
    
    if individual_thread_data:
        print("INDIVIDUAL THREAD USAGE:")
        total_threads = sum(len(threads) for threads in individual_thread_data.values())
        all_times = []
        for rank, threads in individual_thread_data.items():
            for thread_id, thread_time in threads.items():
                thread_time_ms = thread_time / 1000000000.0
                all_times.append(thread_time_ms)
                print(f"  Rank {rank} Thread {thread_id}: {thread_time_ms:.1f} s")
        
        if all_times:
            total_thread_time = sum(all_times)
            max_thread_time = max(all_times)
            min_thread_time = min(all_times)
            balance_ratio = min_thread_time / max_thread_time if max_thread_time > 0 else 0
            
            print(f"  Total threads active: {total_threads}")
            print(f"  Overall total thread time: {total_thread_time:.1f} s")
            print(f"  Thread balance ratio: {balance_ratio:.2f}")
    else:
        print("INDIVIDUAL THREAD USAGE: No data available")
    
    print()
    print(f"Output saved to: pa_vis.pdf")
    print("="*60)

if __name__ == "__main__":
    main()