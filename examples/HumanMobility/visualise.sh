#!/bin/bash

rm images-rivers-3/* # don't remove the images directory if locally

num_files=18
min_x=4000
max_x=6000
min_y=4000
max_y=6000
type="gas"  # or "particles"
river_base="river-rivers-3"
data_base="data-rivers-3/humanMobility"
image_base="images-rivers-3/humanMobility"

# Plot the result
python3 plot_velocity_parallel.py ${num_files} ${min_x} ${max_x} ${min_y} ${max_y} ${type} ${river_base} ${data_base} ${image_base}
# python3 generate_GIF.py ${num_files}

# This command sets:
# - A framerate of 20 frames per second.
# - The input pattern to 'humanMobility_%04d.png' starting from 'humanMobility_0000.png'.
# - The total number of frames to 75 (covering humanMobility_0000.png to humanMobility_0074.png).
# - The output video format to H.264 with yuv420p pixel format for broader compatibility.

# ffmpeg -framerate 100 \ # uncomment this line if locally
#        -start_number 0 \
#        -i cosma/images/humanMobility_%04d.png \
#        -frames:v ${num_files} \
#        -y \
#        -c:v libx264 \
#        -pix_fmt yuv420p \
#        cosma/video.mp4
