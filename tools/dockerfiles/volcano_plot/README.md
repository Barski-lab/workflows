# Interactive Visual Plugins
Interactive visual plugins is a project turned into a docker image to be used in [workflows](https://github.com/datirium/workflows)

The interactive html page is currently only compatible with the volcano_plot, but other plot types will eventually be integrated as well.

- [Purpose](#purpose)
- [How scidap uses plugins](#how-scidap-uses-interactive-plugins)
- [How to customize plugin for scidaps worflows](#how-to-customize-interactive-plugins-for-scidaps-workflows)
- [Example implementation in custom workflow](#custom-workflows-example)


## Development

To view the html report and test it's interactivitiy with js, the html must be served.
Run a python simple http server from within the html_data folder

```bash
# from root dir
python -m http.server 8000 --bind 127.0.0.1
# access at https://localhost:8000/PLOT_FOLDER_NAME/html_data/index.html
```

---

## Purpose

The purpose of this project/repo is to formalize what plots will be generated in the html page, as well as what sample data will be used.
Another goal is to make visual plugins more customizable. Not only through scidap by users (to customize plots), but in developement by devs/BI's (to customize the interactive page itself)

---

## How SciDAP Uses Interactive Plugins

Usage of interactive pages within Scidap will be facilitated by custom visual plugins. These plugins are made available as output files on workflows. 

The logical flows is as follows:
1. workflows run sample analysis in steps
2. one of these steps will generate html (using the docker image built from this project as the template)
3. this step will include a [script](/setVars.sh) that will set the html template variables according to the specific sample/workflow running 
4. viewing the output file through scidap will open an interactive page. the plot type (currently only volcano supported) and format of the page are determined by the html generated on some workflow step



## How to Customize Interactive Plugins for SciDAPs workflows
There are 2 ways to customize the interactive plot. Through [formatting](#format-plot), and through [data setting](#data-setting)

- Customizing format means changing how the form and plot are visually oriented and sized on the page
- Customizing data means changing what file and columns are used for the plot

### Format Plot
Visual plugins are generated on the html template based on div ID's in the template
The div ID of ```sidebar``` is where the form will be created.
The div ID of ```chart``` is where the plot itself will be created

The default format of the page has the form on the left of the page (as a sidebar) and the chart as the content occupying the space remaining. They are formatted using flex layouts, so smaller screens will still have a majority of the screen given to the plot, with a side portion given to the form.

Customizing the layout of the html page requires utilizing [custom workflows in scidap](TODO:add-contributing-how-to-to-scidap)

### Customize Data Used
Deciding what output file will be used (as well as data columns) is established within the workflow step that creates the html page.
Using SciDAP's html template, a script will set variables in the html according to the sample/workflow it is generated for

// TODO: include examples of how to give script: names of data columns, name of file

#### If customizing data in a custom html template
In order for our script to properly set the html template according to the sample/workflow, your custom html template should include (ONLY) the following script in it's < body > tag

```html
<!-- ... -->
<script type="text/javascript">
    var file = window.location.href.replace('index.html', '') + '/report.tsv'; 
    var xColName = 'log2FoldChange';
    var yColName = 'padj';
    var dataNameCol = 'GeneId';
    let plotDivID = "chart";
    let formDivID = "sidebar"; 

    initVisualPlugin(file, plotDivID, formDivID, xColName, yColName, dataNameCol);
</script>
```

NOTE: the variables default values for these variables ARE important, as they are how our script knows what lines and where in the line to include sample/workflow specific data.
---

## Custom Workflows Example
In custom workflows, any html file can be used as an output in order to utilize the interactive visual plugins.
As an output, we want some other steps output to be the data used, so the step creating your custom html page should be downstream to the step creating the data.

// TODO: examples of cwl outputs and html files that can utilize our interactive plugins
// TODO: formalize what the 'sd:visualPlugins' type will be for interactive pages (different type for each plot?, attribute for handling different plots?)