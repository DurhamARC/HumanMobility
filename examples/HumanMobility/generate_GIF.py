import imageio.v2 as imageio

images = []
for i in range(450):
    filename = "images/humanMobility_%04d.png" % i
    images.append(imageio.imread(filename))

# Save as GIF
imageio.mimsave('animation.gif', images, duration=0.1)
