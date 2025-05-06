import numpy as np
import h5py as h5
import argparse as ap
from scipy.ndimage import gaussian_filter
from scipy.interpolate import splprep, splev

def generate_meandering_river(box_size, num_control_points=8, river_width=60.0, randomness=0.15):
    """Generate a meandering river using control points and spline interpolation"""
    
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
    box_size = [args.box_size, args.box_size]  # square domain
    mass = 1e4                    # "mass" of the river
    distance = 1.0                # minimal distance from river
    river_width = 200.0           # width of the river
    
    # Grid parameters
    grid_size = [args.grid_size+1, args.grid_size+1]  # square grid
    
    # Generate meandering river
    x_center, y_center, left_bank_x, left_bank_y, right_bank_x, right_bank_y = \
        generate_meandering_river(box_size, river_width=river_width)
    
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
        f["Header"].attrs["RandomSeed"] = np.array(args.seed, dtype=np.int32)

if __name__ == "__main__":
    main()