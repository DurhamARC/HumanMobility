#!/bin/bash

# Remove previously generated river images (do not remove the directory itself)
rm images-rivers-3/rivers.png 2>/dev/null

# Set parameters
min_x=0
max_x=100000
min_y=0
max_y=100000
river_base="river-rivers-3"
output_image_base="images-rivers-3/rivers"

# Plot only the rivers
python3 plot_rivers.py ${river_base}.hdf5 ${output_image_base}.png ${min_x} ${max_x} ${min_y} ${max_y}