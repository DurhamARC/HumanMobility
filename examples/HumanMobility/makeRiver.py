import numpy as np
import h5py as h5
import argparse as ap
from scipy.ndimage import gaussian_filter

def generate_river(box_size, river_width=60.0):
    """Generate a horizontal river with constant width"""
    
    # Generate river centerline (straight line from west to east)
    x_center = np.linspace(0, box_size[0], 1000)
    y_center = np.ones_like(x_center) * box_size[1]/2
    
    # Generate river banks by offsetting from centerline
    half_width = river_width / 2.0
    left_bank_x = x_center
    left_bank_y = y_center - half_width
    right_bank_x = x_center
    right_bank_y = y_center + half_width
    
    return x_center, y_center, left_bank_x, left_bank_y, right_bank_x, right_bank_y

def calculate_river_acceleration(x, y, left_bank_x, left_bank_y, right_bank_x, right_bank_y, mass, distance):
    """Calculate acceleration at point (x,y) due to river banks
    
    For horizontal river:
    - Inside river: acceleration proportional to distance from centerline
    - Within 'distance' from banks: acceleration proportional to distance constant
    - Outside: acceleration proportional to actual distance from closest bank
    """
    
    # For horizontal river, we only need y-distances
    dy_left = y - left_bank_y[0]  # y-coordinate is constant for each bank
    dy_right = y - right_bank_y[0]
    river_center = (left_bank_y[0] + right_bank_y[0]) / 2.0
    
    # Calculate acceleration
    ax = 0.0
    ay = 0.0
    
    # Inside river check
    if dy_left >= 0 and dy_right <= 0:
        # Inside river - acceleration proportional to distance from centerline
        dy_center = y - river_center
        r = distance  # Use constant distance for magnitude
        rinv3 = 1.0 / (r * r * r)
        # Direction away from centerline
        ay = mass * np.sign(dy_center) * distance * rinv3
    else:
        # Outside river or near banks
        # Find distance to closest bank
        if abs(dy_left) < abs(dy_right):
            # Closer to left bank
            if abs(dy_left) <= distance:
                # Within distance constant from bank - use constant force
                r = distance
                rinv3 = 1.0 / (r * r * r)
                ay = mass * np.sign(dy_left) * distance * rinv3
            else:
                # Beyond distance constant - use actual distance
                r = abs(dy_left)
                rinv3 = 1.0 / (r * r * r)
                ay = mass * dy_left * rinv3
        else:
            # Closer to right bank
            if abs(dy_right) <= distance:
                # Within distance constant from bank - use constant force
                r = distance
                rinv3 = 1.0 / (r * r * r)
                ay = mass * np.sign(dy_right) * distance * rinv3
            else:
                # Beyond distance constant - use actual distance
                r = abs(dy_right)
                rinv3 = 1.0 / (r * r * r)
                ay = mass * dy_right * rinv3
    
    return ax, ay

def main():
    # Parse arguments
    parser = ap.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, default="river.hdf5")
    parser.add_argument("-s", "--seed", type=int, default=42,
                       help="Random seed for river generation")
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    np.random.seed(args.seed)
    
    # Parameters
    box_size = [10000., 10000.]  # meters
    mass = 1e5                    # "mass" of the river
    distance = 1.0                # minimal distance from river
    river_width = 60.0           # width of the river
    
    # Grid parameters
    grid_size = [1001, 1001]     # number of cells + 1
    
    # Generate straight river
    x_center, y_center, left_bank_x, left_bank_y, right_bank_x, right_bank_y = \
        generate_river(box_size, river_width=river_width)
    
    # Create acceleration field arrays
    ax = np.zeros(grid_size, dtype=np.float32)
    ay = np.zeros(grid_size, dtype=np.float32)
    
    # Generate grid coordinates
    x_coords = np.linspace(0, box_size[0], grid_size[0])
    y_coords = np.linspace(0, box_size[1], grid_size[1])
    
    # Calculate acceleration field at each grid point
    print("Calculating acceleration field...")
    for i, x in enumerate(x_coords):
        if i % 100 == 0:  # Progress indicator
            print(f"Processing row {i}/{grid_size[0]}")
        for j, y in enumerate(y_coords):
            ax[i,j], ay[i,j] = calculate_river_acceleration(
                x, y, left_bank_x, left_bank_y, right_bank_x, right_bank_y,
                mass, distance
            )
    
    # Optional: Smooth the acceleration field
    # print("Smoothing acceleration field...")
    # ax = gaussian_filter(ax, sigma=1.0)
    # ay = gaussian_filter(ay, sigma=1.0)
    
    # Save to HDF5 file
    print("Saving to HDF5 file...")
    with h5.File(args.file, 'w') as f:
        # Save acceleration field
        f.create_dataset("AccelerationField/ax", data=ax)
        f.create_dataset("AccelerationField/ay", data=ay)
        
        # Save river geometry for visualization
        f.create_dataset("RiverGeometry/centerline_x", data=x_center)
        f.create_dataset("RiverGeometry/centerline_y", data=y_center)
        f.create_dataset("RiverGeometry/left_bank_x", data=left_bank_x)
        f.create_dataset("RiverGeometry/left_bank_y", data=left_bank_y)
        f.create_dataset("RiverGeometry/right_bank_x", data=right_bank_x)
        f.create_dataset("RiverGeometry/right_bank_y", data=right_bank_y)
        
        # Save metadata
        f.create_group("Header")
        f["Header"].attrs["BoxSize"] = np.array(box_size, dtype=np.float64)
        f["Header"].attrs["GridSize"] = np.array(grid_size, dtype=np.int32)
        f["Header"].attrs["Mass"] = np.array(mass, dtype=np.float64)
        f["Header"].attrs["Distance"] = np.array(distance, dtype=np.float64)
        f["Header"].attrs["RiverWidth"] = np.array(river_width, dtype=np.float64)
        f["Header"].attrs["RandomSeed"] = np.array(args.seed, dtype=np.int32)  # Add this line

if __name__ == "__main__":
    main()
