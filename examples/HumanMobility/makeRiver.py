import numpy as np
import h5py as h5
import argparse as ap
from scipy.ndimage import gaussian_filter

# Parse arguments
parser = ap.ArgumentParser()
parser.add_argument("-f", "--file", type=str, default="river.hdf5")
args = parser.parse_args()

# River parameters (matching humanMobility.yml)
box_size = [10000., 10000.]  # meters
river_y = [5010., 5090.]     # river banks y-coordinates
river_width = river_y[1] - river_y[0]
river_center = np.mean(river_y)

# Grid parameters
grid_size = [100, 100]  # number of cells
dx = box_size[0] / grid_size[0]
dy = box_size[1] / grid_size[1]

# Create acceleration field
ax = np.zeros(grid_size)
ay = np.zeros(grid_size)

# Generate grid coordinates
y_coords = np.linspace(0, box_size[1], grid_size[1])

# Set acceleration field
for i in range(grid_size[1]):
    y = y_coords[i]
    # Inside river
    if river_y[0] <= y <= river_y[1]:
        ax[:, i] = 1.0  # Flow along x-direction
    else:
        # Decay outside river
        dist = min(abs(y - river_y[0]), abs(y - river_y[1]))
        decay = np.exp(-dist / river_width)
        ax[:, i] = decay

# Smooth the acceleration field
ax = gaussian_filter(ax, sigma=2.0)
ay = gaussian_filter(ay, sigma=2.0)

# Save to HDF5 file
with h5.File(args.file, 'w') as f:
    f.create_dataset("AccelerationField/ax", data=ax)
    f.create_dataset("AccelerationField/ay", data=ay)
    f.create_dataset("Header/BoxSize", data=box_size)
    f.create_dataset("Header/GridSize", data=grid_size)