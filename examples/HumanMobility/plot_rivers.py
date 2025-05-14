import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import os

plt.style.use("../../tools/stylesheets/mnras.mplstyle")

def plot_rivers(river_filename, output_png, min_x=None, max_x=None, min_y=None, max_y=None):
    with h5py.File(river_filename, 'r') as sim:
        box_size = sim["/Header"].attrs["BoxSize"]
        river_group = sim["/RiverGeometry"]
        num_rivers = river_group.attrs.get("num_rivers", 1)
        river_geometries = []
        for idx in range(num_rivers):
            grp = river_group[f"river_{idx}"]
            river_geometries.append({
                "left_bank_x": grp["left_bank_x"][:],
                "left_bank_y": grp["left_bank_y"][:],
                "right_bank_x": grp["right_bank_x"][:],
                "right_bank_y": grp["right_bank_y"][:],
                "center_x": grp["centerline_x"][:],
                "center_y": grp["centerline_y"][:],
            })

    fig, ax = plt.subplots()
    ax.set_title("River Geometry")

    # Plot river banks and centerline
    for river in river_geometries:
        ax.plot(river["left_bank_x"], river["left_bank_y"], '-', color='b', linewidth=0.8, label=None)
        ax.plot(river["right_bank_x"], river["right_bank_y"], '-', color='b', linewidth=0.8, label=None)
        ax.plot(river["center_x"], river["center_y"], '--', color='k', linewidth=0.5, label=None)

    if min_x is not None and max_x is not None:
        ax.set_xlim(min_x, max_x)
    elif "BoxSize" in locals():
        ax.set_xlim(0, box_size[0])
    if min_y is not None and max_y is not None:
        ax.set_ylim(min_y, max_y)
    elif "BoxSize" in locals():
        ax.set_ylim(0, box_size[1])

    plt.xlabel("${\\rm{Position}}~x$")
    plt.ylabel("${\\rm{Position}}~y$")
    plt.savefig(output_png, dpi=300)
    plt.close()

if __name__ == '__main__':
    # Usage: python plot_rivers_only.py river_file.hdf5 output.png [min_x max_x min_y max_y]
    if len(sys.argv) < 3:
        print("Usage: python plot_rivers_only.py river_file.hdf5 output.png [min_x max_x min_y max_y]")
        sys.exit(1)
    river_filename = sys.argv[1]
    output_png = sys.argv[2]
    min_x = float(sys.argv[3]) if len(sys.argv) > 3 else None
    max_x = float(sys.argv[4]) if len(sys.argv) > 4 else None
    min_y = float(sys.argv[5]) if len(sys.argv) > 5 else None
    max_y = float(sys.argv[6]) if len(sys.argv) > 6 else None

    plot_rivers(river_filename, output_png, min_x, max_x, min_y, max_y)