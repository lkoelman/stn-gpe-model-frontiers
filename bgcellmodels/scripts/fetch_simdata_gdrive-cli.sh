#!/bin/bash -l

# Using google drive client github.com/odeke-em/drive

# First push new data to google drive on cluster
# > drive push

# Filter data using regular expression and print separated by spaces
# Also remove the leading forward slash on each line
pulldata=$(drive ls simdata_newsonic | awk 'BEGIN { ORS=" " } /2019.08.01/ { print substr($0,2) }')

# Pull matching data
drive pull $pulldata


# Alternative: put data in new subfolder and pull only that
drive pull -no-clobber simdata_newsonic/subfolder

