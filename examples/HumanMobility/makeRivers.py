import numpy as np
import h5py as h5
import argparse as ap
from scipy.ndimage import gaussian_filter
from scipy.interpolate import splprep, splev
import concurrent.futures
import time
import os

def generate_meandering_river(box_size, num_control_points=8, river_width=60.0, randomness=0.15):
    """Generate a meandering river using control points and spline interpolation"""

    print("Generating meandering river...")
    
    # Generate random control points for river centerline
    # Start from left side (west)
    x_controls = np.linspace(0, box_size[0], num_control_points)
    
    # Add random vertical displacement to control points (keeping ends fixed)
    y_controls = np.zeros(num_control_points)
    y_controls[1:-1] = box_size[1]/2 + randomness * box_size[1] * np.random.randn(num_control_points-2)
    y_controls[0] = box_size[1]/2  # Fix start point
    y_controls[-1] = box_size[1]/2  # Fix end point
    
    # Fit a spline through control points
    tck, u = splprep([x_controls, y_controls], s=0, k=3)
    
    # Generate points along the spline for smooth river
    t_fine = np.linspace(0, 1, 1000)
    x_center, y_center = splev(t_fine, tck)
    
    # Calculate tangent vectors along the river
    dx = np.gradient(x_center)
    dy = np.gradient(y_center)
    
    # Normalize tangent vectors
    norm = np.sqrt(dx*dx + dy*dy)
    dx /= norm
    dy /= norm
    
    # Calculate normal vectors (perpendicular to tangent)
    normal_x = -dy
    normal_y = dx
    
    # Generate river banks by offsetting perpendicular to centerline
    half_width = river_width / 2.0
    left_bank_x = x_center - normal_x * half_width
    left_bank_y = y_center - normal_y * half_width
    right_bank_x = x_center + normal_x * half_width
    right_bank_y = y_center + normal_y * half_width
    
    return x_center, y_center, left_bank_x, left_bank_y, right_bank_x, right_bank_y

def calculate_river_acceleration(x, y, left_bank_x, left_bank_y, right_bank_x, right_bank_y, mass, distance):
    """Calculate acceleration at point (x,y) due to river banks
    
    For meandering river:
    - Inside river: acceleration proportional to distance from centerline
    - Within 'distance' from banks: acceleration proportional to distance constant
    - Outside: acceleration proportional to actual distance from closest bank
    """

    # Find distances and closest points on both banks
    left_dists = np.sqrt((x - left_bank_x)**2 + (y - left_bank_y)**2)
    right_dists = np.sqrt((x - right_bank_x)**2 + (y - right_bank_y)**2)
    
    left_idx = np.argmin(left_dists)
    right_idx = np.argmin(right_dists)
    min_left_dist = left_dists[left_idx]
    min_right_dist = right_dists[right_idx]
    
    # Calculate vectors from banks to point
    dx_left = x - left_bank_x[left_idx]
    dy_left = y - left_bank_y[left_idx]
    dx_right = x - right_bank_x[right_idx]
    dy_right = y - right_bank_y[right_idx]

    # Calculate normal vectors to check if point is between banks
    normal_x = right_bank_x[right_idx] - left_bank_x[left_idx]
    normal_y = right_bank_y[right_idx] - left_bank_y[left_idx]
    norm = np.sqrt(normal_x*normal_x + normal_y*normal_y)
    normal_x /= norm
    normal_y /= norm

    ax = 0.0
    ay = 0.0

    # Inside river check (if between banks)
    dot_left = dx_left*normal_x + dy_left*normal_y
    dot_right = dx_right*normal_x + dy_right*normal_y
    if dot_left >= 0 and dot_right <= 0:  # This checks if point is between banks
        # Inside river - acceleration proportional to distance from centerline
        center_x = (left_bank_x[left_idx] + right_bank_x[right_idx]) / 2
        center_y = (left_bank_y[left_idx] + right_bank_y[right_idx]) / 2
        dx = x - center_x
        dy = y - center_y
        r = distance  # Use constant distance for magnitude
        rinv3 = 1.0 / (r * r * r)
        # Direction away from centerline
        norm = np.sqrt(dx*dx + dy*dy)
        if norm > 0:
            ax = mass * (dx/norm) * distance * rinv3
            ay = mass * (dy/norm) * distance * rinv3
    else:
        # Outside river or near banks
        # Find closest bank
        if min_left_dist < min_right_dist:
            # Closer to left bank
            dx = x - left_bank_x[left_idx]
            dy = y - left_bank_y[left_idx]
            if min_left_dist <= distance:
                # Within distance constant from bank - use constant force
                r = distance
                rinv3 = 1.0 / (r * r * r)
                norm = np.sqrt(dx*dx + dy*dy)
                ax = mass * (dx/norm) * distance * rinv3
                ay = mass * (dy/norm) * distance * rinv3
            else:
                # Beyond distance constant - use actual distance
                r = min_left_dist
                rinv3 = 1.0 / (r * r * r)
                norm = np.sqrt(dx*dx + dy*dy)
                ax = mass * (dx/norm) * r * rinv3
                ay = mass * (dy/norm) * r * rinv3
        else:
            # Closer to right bank
            dx = x - right_bank_x[right_idx]
            dy = y - right_bank_y[right_idx]
            if min_right_dist <= distance:
                # Within distance constant from bank - use constant force
                r = distance
                rinv3 = 1.0 / (r * r * r)
                norm = np.sqrt(dx*dx + dy*dy)
                ax = mass * (dx/norm) * distance * rinv3
                ay = mass * (dy/norm) * distance * rinv3
            else:
                # Beyond distance constant - use actual distance
                r = min_right_dist
                rinv3 = 1.0 / (r * r * r)
                norm = np.sqrt(dx*dx + dy*dy)
                ax = mass * (dx/norm) * r * rinv3
                ay = mass * (dy/norm) * r * rinv3
    
    return ax, ay

def compute_row_block(i, x, y_coords, river_segments, mass, distance):
    ax_row = np.zeros(len(y_coords), dtype=np.float32)
    ay_row = np.zeros(len(y_coords), dtype=np.float32)
    for j, y in enumerate(y_coords):
        min_dist = float('inf')
        best_seg = None
        for seg in river_segments:
            seg_y = np.mean(seg["centerline_y"])
            seg_x = np.mean(seg["centerline_x"])
            dist = np.sqrt((x - seg_x)**2 + (y - seg_y)**2)
            if dist < min_dist:
                min_dist = dist
                best_seg = seg
        ax_row[j], ay_row[j] = calculate_river_acceleration(
            x, y,
            best_seg["left_bank_x"], best_seg["left_bank_y"],
            best_seg["right_bank_x"], best_seg["right_bank_y"],
            mass, distance
        )
    return i, ax_row, ay_row

def main():
    # Parse arguments
    parser = ap.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, default="river.hdf5",
                       help="Output HDF5 file name")
    parser.add_argument("-s", "--seed", type=int, default=42,
                       help="Random seed for river generation")
    parser.add_argument("-b", "--box-size", type=float, default=10000.,
                       help="Box size in meters (square domain)")
    parser.add_argument("-g", "--grid-size", type=int, default=1000,
                       help="Number of grid cells in each dimension + 1")
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    np.random.seed(args.seed)
    
    # Parameters
    max_box = min(args.box_size, 10000.0)
    box_size = [max_box, max_box]  # restrict to [0, min(args.box_size, 10000)]
    mass = 1e4                    # "mass" of the river
    distance = 1.0                # minimal distance from river
    river_width = 200.0           # width of the river
    
    # Grid parameters
    grid_size = [args.grid_size+1, args.grid_size+1]  # square grid

    river_segments = []
    block_size = 10000.0
    start_time = time.time()
    if args.box_size > block_size:
        n_repeat_x = int(np.ceil(args.box_size / block_size))
        n_repeat_y = int(np.ceil(args.box_size / block_size))

        # Generate the base river only once in the bottom left corner
        base_box_size = [block_size, block_size]
        x_c, y_c, l_x, l_y, r_x, r_y = generate_meandering_river(base_box_size, river_width=river_width)

        # Replicate the river along X and Y axes
        for j in range(n_repeat_y):
            for i in range(n_repeat_x):
                offset_x = i * block_size
                offset_y = j * block_size
                # Only keep points within the block (should always be true for the base river)
                mask = (x_c >= 0) & (x_c <= block_size)
                river_segments.append({
                    "centerline_x": x_c[mask] + offset_x,
                    "centerline_y": y_c[mask] + offset_y,
                    "left_bank_x": l_x[mask] + offset_x,
                    "left_bank_y": l_y[mask] + offset_y,
                    "right_bank_x": r_x[mask] + offset_x,
                    "right_bank_y": r_y[mask] + offset_y,
                })
        box_size = [args.box_size, args.box_size]
    else:
        x_center, y_center, left_bank_x, left_bank_y, right_bank_x, right_bank_y = \
            generate_meandering_river(box_size, river_width=river_width)
        river_segments.append({
            "centerline_x": x_center,
            "centerline_y": y_center,
            "left_bank_x": left_bank_x,
            "left_bank_y": left_bank_y,
            "right_bank_x": right_bank_x,
            "right_bank_y": right_bank_y,
        })
    end_time = time.time()
    print(f"River generation took {end_time - start_time:.2f} seconds.")

    # Extended information output
    print(f"Block size: {block_size} x {block_size}")
    print(f"Grid size in blocks: {n_repeat_x} x {n_repeat_y}" if args.box_size > block_size else "Grid size in blocks: 1 x 1")
    print(f"Full domain size: {box_size[0]} x {box_size[1]} (meters)")
    print(f"Full grid size: {grid_size[0]} x {grid_size[1]}")

    # Calculate number of blocks
    n_repeat_x = int(np.ceil(args.box_size / block_size))
    n_repeat_y = int(np.ceil(args.box_size / block_size))

    # Grid parameters
    grid_size = [args.grid_size+1, args.grid_size+1]  # e.g. 10001 x 10001

    # Calculate block grid size (number of points per block, including overlap)
    block_grid_size_x = (grid_size[0] - 1) // n_repeat_x + 1
    block_grid_size_y = (grid_size[1] - 1) // n_repeat_y + 1
    block_grid_size = [block_grid_size_x, block_grid_size_y]

    # Generate grid coordinates for the block
    x_block = np.linspace(0, block_size, block_grid_size[0])
    y_block = np.linspace(0, block_size, block_grid_size[1])

    # Create acceleration field arrays for a single block
    ax_block = np.zeros(block_grid_size, dtype=np.float32)
    ay_block = np.zeros(block_grid_size, dtype=np.float32)

    # Calculate acceleration field for the bottom-left block
    print("Calculating acceleration field for bottom-left block...")
    accel_start_time = time.time()

    # Limit the number of worker processes to avoid "too many open files" error
    # Use at most 16 processes or the number of CPU cores, whichever is smaller
    max_workers = min(16, os.cpu_count() or 1)
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for i, x in enumerate(x_block):
            if(i % 100 == 0):
                print(f"Processing row {i+1}/{block_grid_size[0]}")
            # Submit tasks to the executor
            futures.append(executor.submit(compute_row_block, i, x, y_block, river_segments, mass, distance))
        
        for future in concurrent.futures.as_completed(futures):
            i, ax_row, ay_row = future.result()
            ax_block[i, :] = ax_row
            ay_block[i, :] = ay_row
    accel_end_time = time.time()
    print(f"Acceleration field calculation for block took {accel_end_time - accel_start_time:.2f} seconds.")

    # Now assemble the full grid, avoiding duplicated edges
    ax = np.zeros(grid_size, dtype=np.float32)
    ay = np.zeros(grid_size, dtype=np.float32)

    for i in range(n_repeat_x):
        for j in range(n_repeat_y):
            # Compute start and end indices for this block in the full grid
            x_start = i * (block_grid_size_x - 1)
            x_end = x_start + block_grid_size_x
            y_start = j * (block_grid_size_y - 1)
            y_end = y_start + block_grid_size_y

            # For the last block, ensure we don't go out of bounds
            if x_end > grid_size[0]:
                x_end = grid_size[0]
            if y_end > grid_size[1]:
                y_end = grid_size[1]

            # Compute the corresponding slice in the block grid
            bx_start = 0
            bx_end = x_end - x_start
            by_start = 0
            by_end = y_end - y_start

            ax[x_start:x_end, y_start:y_end] = ax_block[bx_start:bx_end, by_start:by_end]
            ay[x_start:x_end, y_start:y_end] = ay_block[bx_start:bx_end, by_start:by_end]

    # Save to HDF5 file
    print("Saving to HDF5 file...")
    with h5.File(args.file, 'w') as f:
        # Save acceleration field
        f.create_dataset("AccelerationField/ax", data=ax)
        f.create_dataset("AccelerationField/ay", data=ay)
        
        # Save river geometry for visualization as a collection
        river_group = f.create_group("RiverGeometry")
        for idx, seg in enumerate(river_segments):
            grp = river_group.create_group(f"river_{idx}")
            grp.create_dataset("centerline_x", data=seg["centerline_x"])
            grp.create_dataset("centerline_y", data=seg["centerline_y"])
            grp.create_dataset("left_bank_x", data=seg["left_bank_x"])
            grp.create_dataset("left_bank_y", data=seg["left_bank_y"])
            grp.create_dataset("right_bank_x", data=seg["right_bank_x"])
            grp.create_dataset("right_bank_y", data=seg["right_bank_y"])
        river_group.attrs["num_rivers"] = len(river_segments)
        
        # Save metadata
        f.create_group("Header")
        f["Header"].attrs["BoxSize"] = np.array(box_size, dtype=np.float64)
        f["Header"].attrs["GridSize"] = np.array(grid_size, dtype=np.int32)
        f["Header"].attrs["Mass"] = np.array(mass, dtype=np.float64)
        f["Header"].attrs["Distance"] = np.array(distance, dtype=np.float64)
        f["Header"].attrs["RiverWidth"] = np.array(river_width, dtype=np.float64)
        f["Header"].attrs["RandomSeed"] = np.array(args.seed, dtype=np.int32)

if __name__ == "__main__":
    main()