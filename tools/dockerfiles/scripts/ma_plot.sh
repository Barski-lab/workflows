#!/bin/bash

echo "Copying volcano plot to the current directory"
cp -r /opt/volcano_plot .
cd ./volcano_plot
./MA_PLOT_set_vars.sh "$(basename -- $1)" $2 $3 $4 "chart" "sideForm" "MA" "MD-MA_plot/html_data"