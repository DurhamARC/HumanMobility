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
mass = 1e6                 # "mass" of the river
distance = 1.0          # minimal distance from the river

# Grid parameters
x_grid_size = 1000
y_grid_size = 1000
grid_size = [x_grid_size+1, y_grid_size+1]  # number of cells

# Create acceleration field
ax = np.zeros(grid_size, dtype=np.float32)
ay = np.zeros(grid_size, dtype=np.float32)

# Generate grid coordinates
y_coords = np.linspace(0, box_size[1], grid_size[1])

# Set acceleration field
# ax component is zero outside river
ax[:, :] = 0
for i in range(grid_size[1]):
    y = y_coords[i]

### FOR THE RIVER FROM WEST TO EAST (simple example) ###
    # Determine distance from river
    if y < river_y[0]-distance:
        dy = y - river_y[0]
    elif y > river_y[1]+distance:
        dy = y - river_y[1]
    else:
        # Inside the river
        dy = distance
    
    # Calculate acceleration using inverse cube law (rinv3)
    r = np.sqrt(dy * dy)
    rinv = 1.0 / r
    rinv3 = rinv * rinv * rinv
    
    # Set acceleration components
    # ay component follows inverse cube law
    ay[:, i] = mass * dy * rinv3

# Smooth the acceleration field
#ax = gaussian_filter(ax, sigma=1.0)
#ay = gaussian_filter(ay, sigma=1.0)
### END RIVER ###

# Save to HDF5 file
with h5.File(args.file, 'w') as f:
    print(f"Creating acceleration field with shape: {ax.shape}")
    print(f"Grid size: {grid_size}")
    print(f"Box size: {box_size}")
    print(f"Mass: {mass}")
    
    f.create_dataset("AccelerationField/ax", data=ax)
    f.create_dataset("AccelerationField/ay", data=ay)
    f.create_group("Header")
    f["Header"].attrs["BoxSize"] = np.array(box_size, dtype=np.float64)
    f["Header"].attrs["GridSize"] = np.array(grid_size, dtype=np.int32)
    f["Header"].attrs["Mass"] = np.array(mass, dtype=np.float64)
    f["Header"].attrs["Distance"] = np.array(distance, dtype=np.float64)
    
    # Verify data was written correctly
    print(f"Wrote ax with shape: {f['AccelerationField/ax'].shape}")
    print(f"Wrote ay with shape: {f['AccelerationField/ay'].shape}")
    print(f"Header/BoxSize: {f['Header'].attrs['BoxSize']}")
    print(f"Header/GridSize: {f['Header'].attrs['GridSize']}")
    print(f"Header/GridSize: {f['Header'].attrs['Mass']}")