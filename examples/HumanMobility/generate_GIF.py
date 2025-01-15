import sys
import imageio.v2 as imageio

if __name__ == '__main__':

    num_files = int(sys.argv[1])   # Adjust to the number of your PNG files
    images = []
    for i in range(num_files):
        filename = "images/humanMobility_%04d.png" % i
        images.append(imageio.imread(filename))

    # Save as GIF
    imageio.mimsave('animation.gif', images, duration=0.1)
