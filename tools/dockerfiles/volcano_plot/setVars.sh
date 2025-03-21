#!/bin/bash

# script to set the 4 main variables used to create the plot
# $1: name of tsv file holding data for the plot (ex: "report.tsv")
# $2: name of column in file for the plots x-axis (ex: "log2FoldChange")
# $3: name of column in file for the plots y-axis (ex: "padj")
# $4: name of column in file for each data points 'name' (ex: "GeneId")
# $5: ID of div in html that will hold the plot (default: 'sidebar')
# $6: ID of div in html that will hold the plot (default: 'chart')
# $7: type of plot to be created (ex: 'volcano')
# $8: name of parent folder this data is stored in on the satellite

# test script with following command (run from repo root)
# ./setVars.sh report_new.tsv new_x_xol new_y_col new_data_col sidebar chart volcano

## default usage of script
## ./setVars.sh volcano_report.tsv log2FoldChange padj GeneId chart sidebar volcano volcano_plot/html_data

DATA_FILE=$1
X_COL_NAME=$2
Y_COL_NAME=$3
DATA_NAME_COL=$4
PLOT_DIV_ID=$5
FORM_DIV_ID=$6
PLOT_TYPE=$7
PARENT_FOLDER=$8
INDEX_HTML_FILE=$9

# for each variable to set according to workflow: create new line for html -> sed command with old line and new line

NEW_DATA_FILE_LINE="var file = window.location.href.replace('\/$PARENT_FOLDER' + '\/$INDEX_HTML_FILE.html', '') + '\/$DATA_FILE';"
#echo $NEW_DATA_FILE_LINE
sed -i "s@var file = window.location.href.replace('index.html', '') + '\/report.tsv';@$NEW_DATA_FILE_LINE@" volcano_plot/html_data/index.html

NEW_X_COL_LINE="var xColName = '$X_COL_NAME';"
#echo $NEW_X_COL_LINE
sed -i "s/var xColName = 'log2FoldChange';/$NEW_X_COL_LINE/" volcano_plot/html_data/index.html

NEW_Y_COL_LINE="var yColName = '$Y_COL_NAME';"
#echo $NEW_Y_COL_LINE
sed -i "s/var yColName = 'padj';/$NEW_Y_COL_LINE/" volcano_plot/html_data/index.html

NEW_DATA_NAME_LINE="var dataNameCol = '$DATA_NAME_COL';"
#echo $NEW_DATA_NAME_LINE
sed -i "s/var dataNameCol = 'GeneId';/$NEW_DATA_NAME_LINE/" volcano_plot/html_data/index.html

NEW_PLOT_DIV_LINE="let plotDivID = '$PLOT_DIV_ID';"
#echo $NEW_PLOT_DIV_LINE
sed -i "s/let plotDivID = 'chart';/$NEW_PLOT_DIV_LINE/" volcano_plot/html_data/index.html

NEW_FORM_DIV_LINE="let formDivID = '$FORM_DIV_ID';"
#echo $NEW_FORM_DIV_LINE
sed -i "s/let formDivID = 'sidebar';/$NEW_FORM_DIV_LINE/" volcano_plot/html_data/index.html

NEW_PLOT_TYPE_LINE="let plotType = '$PLOT_TYPE';"
#echo $NEW_PLOT_TYPE_LINE
sed -i "s/let plotType = 'volcano';/$NEW_PLOT_TYPE_LINE/" volcano_plot/html_data/index.html
