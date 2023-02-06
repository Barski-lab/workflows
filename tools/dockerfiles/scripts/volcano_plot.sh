#!/bin/bash

echo "Copying volcano plot to the current directory"
cp -r /opt/volcano_plot .
cd ./volcano_plot
./setVars.sh "$(basename -- $1)" $2 $3 $4 "chart" "sidebar" "volcano" "volcano_plot/html_data"