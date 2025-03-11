/**
 * 3 functions of importance
 * @method volcanoPlot(): generates volcano plot
 * @method createForm(): creates form for plot
 * @method initVisualPlugins(): adds form and plot to html based on given params
 * 
 * helper methods: parsers for getting data from form values
 */


/**
 * returns a chart that is updated interactively. 
 * data file for chart is set by @method initVisualPlugin() 
 */
 function volcanoPlot() {
    /** below are the variables used to track positioning of different chart components:
     

        width = 960,                                                                        // default for sizing
        height = 500,                                                                       // default for sizing
        margin = { top: 20, right: 20, bottom: 40, left: 50 },                              // default for sizing
        
        xColumn = "log2FoldChange",                                                         // TODO: name column in data file for plotting X-axis
        yColumn = "padj",                                                                   // TODO: name column in data file for plotting Y-axis
        
        xAxisLabel = d3.select("#plotxaxis").property("value"),                             // TODO: label for x axis                                                                                       
        yAxisLabel = d3.select("#plotyaxis").property("value"),                             // TODO: label for y axis
                                                                                            // x&y axis label currently have defaults established in html. moving to here requires some more work

        xAxisLabelOffset,                                                                   // offset for the label of the axis (not customizable)
        yAxisLabelOffset,
        
        xTicks,                                                                             // number of ticks on the axis (not customizable)
        yTicks,
        
        sampleID = "GeneId",                                                                  // TODO: named according to workflow. is column name for getting 'name' of each data point (which gene)       
        
        significanceThreshold = parseFloat(d3.select("#fdrthreshold").property("value")),   // significance threshold to colour by
        foldChangeThreshold = parseFloat(d3.select("#foldchange").property("value")),       // fold change level to colour by
                                                                                            // threshholds currently have default coming from html. more work to have them settable here
        
                                                                                            colorRange,                                                                         // colour range to use in the plot (probably deprecated)

        xScale = d3.scaleLinear(),                                                          // the values for the x-axes will be continuous
        yScale = d3.scaleLog(),                                                             // values for y-axis will be logarithmic
                                                                                            // both should be settable 
        plotTitle = d3.select('#plottitle').property("value"),

                                                                                            // plot dot styles 
        plotDotSize = parseInt(d3.select("#plotdotsize").property("value")), 
        plotDotTransparency = parseFloat(d3.select("#plotdotopacity").property("value")),
        acceptedColorXpos = d3.select("#acceptedxcolorpos").property("value"),
        acceptedColorXneg = d3.select("#acceptedxcolorneg").property("value"),
        rejectedColorX = d3.select("#rejectedxcolor").property("value"),
        rejectedColorY = d3.select("#rejectedycolor").property("value"),
        rejectedColorBoth = d3.select("#rejectedbothcolor").property("value");
        highlightPlotDotSize = parseInt(d3.select("#highlightedplotdotsize").property("value")),
        highlightPlotDotTransparency = parseFloat(d3.select("#highlightedplotdotopacity").property("value")),

        genesToHighlightList = []                                                           // list of data points to higlight in plot
        genesList = []                                                                      // list of all data point names (sampleID col in data)
    
     */
    
    // pretty much every var with a d3.select needs to have work done to make it settable according to workflow
    var //width = "100%",//document.getElementById("temp1").clientWidth, // parent div ID should be given to this function somehow
        //height = "100%",//document.getElementById("temp1").clientHeight,
        width = 960,       // keep static                                                            
        height = 500,                    // change based on custom aspect ratio                                            
        margin = { top: 20, right: 20, bottom: 40, left: 50 },                          
        xColumn = "",                                                     
        yColumn = "",                                                               
        xAxisLabel = d3.select("#plotxaxis").property("value"),
        yAxisLabel = d3.select("#plotyaxis").property("value"),
        xAxisLabelOffset, 
        yAxisLabelOffset,
        xTicks, 
        yTicks,
        sampleID = "", 
        significanceThreshold = parseFloat(d3.select("#fdrthreshold").property("value")), 
        foldChangeThreshold = parseFloat(d3.select("#foldchange").property("value")), 
        colorRange, 
        xScale = d3.scaleLinear(),
        yScale = d3.scaleLog(),
        plotTitle = d3.select('#plottitle').property("value"),
        plotDotSize = parseInt(d3.select("#plotdotsize").property("value")),
        plotDotTransparency = parseFloat(d3.select("#plotdotopacity").property("value")),
        acceptedColorXpos = d3.select("#acceptedxcolorpos").property("value"),
        acceptedColorXneg = d3.select("#acceptedxcolorneg").property("value"),
        rejectedColorX = d3.select("#rejectedxcolor").property("value"),
        rejectedColorY = d3.select("#rejectedycolor").property("value"),
        rejectedColorBoth = d3.select("#rejectedbothcolor").property("value"),
        highlightColor = d3.select("#highlightcolor").property("value"),
        highlightPlotDotSize = parseInt(d3.select("#highlightedplotdotsize").property("value")),
        highlightPlotDotTransparency = parseFloat(d3.select("#highlightedplotdotopacity").property("value")),
        highlightTopXBy = d3.select('#highlighttopn-selections-by').property('value'),
        highlightTopXFilter = d3.select('#highlighttopn-selections-filter').property('value'),
        highlightTopXGenes = parseFloat(d3.select('#highlighttopn').property("value"));//,
        //genesToHighlightList = [];
        //genesList = [];



    function chart(selection) {
        var innerWidth = width - margin.left - margin.right, 
            innerHeight = height - margin.top - margin.bottom; // set the size of the chart according to its container

        selection.each(function (data) {

            // set up the scaling for the axes based on the inner width/height of the chart and also the range
            // of value for the x and y axis variables. This range is defined by their min and max values as
            // calculated by d3.extent()
            xScale.range([0, innerWidth])
                .domain(d3.extent(data, function (d) { return d[xColumn]; }))
                .nice();

            // normally would set the y-range to [height, 0] but by swapping it I can flip the axis and thus
            // have -log10 scale without having to do extra parsing
            yScale.range([0, innerHeight])
                .domain(d3.extent(data, function (d) { return d[yColumn]; }))
                .nice(); // adds "padding" so the domain extent is exactly the min and max values

            /*
                        var zoom = d3.zoom()
                            .scaleExtent([1, 20])
                            .translateExtent([[0, 0], [width, height]])
                            .on('zoom', zoomFunction);
            */
            // append the svg object to the selection. svg object is what is exported as svg/png/pdf
            var svg = d3.select(this).append('svg')
                .attr('id', 'svg')
                //.attr('height', height)
                //.attr('width', width)
                // make height and width come from container but enforce aspect ratio
                .attr("preserveAspectRatio", "xMaxYMax slice")
                .attr("viewBox", `0 0 ${width} ${height}`)   
                //.style('color', 'white') // make background setting part of the 'advanced plot options' tab (not working?)
                .style('background-color', 'white')
                .append('g')
                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');
            //                .call(zoom);

            // position the reset button and attach reset function
            /*d3.select('#resetBtn')
                .style('top', margin.top * 1.5 + 'px')
                .style('left', margin.left * 1.25 + 'px')
                .on('click', resetZoom);*/

            svg.append('defs').append('clipPath')
                .attr('id', 'clip')
                .append('rect')
                .attr('height', innerHeight)
                .attr('width', innerWidth);

            var curPlotTitle = svg.append("text")
                .attr('transform', 'translate(' + width / 2 + ',' + (margin.top - 25) + ')')
                .attr("text-anchor", "middle")  
                .style("font-size", "16px") 
                .style("text-decoration", "underline")  
                .text(plotTitle);

            // add the axes
            var xAxis = d3.axisBottom(xScale);
            var yAxis = d3.axisLeft(yScale)
                .ticks(5)
                .tickFormat(yTickFormat);

            var gX = svg.append('g')
                .attr('class', 'x axis')
                .attr('transform', 'translate(0,' + innerHeight + ')')
                .call(xAxis);

            var gxLabel = gX.append('text')
                .style('fill', '#666')
                .style('font-weight', '700')
                .style('font-size', '12px')
                .attr('transform', 'translate(' + width / 2 + ',' + (margin.bottom - 6) + ')')
                .attr('text-anchor', 'middle')
                .html(xAxisLabel || xColumn);

            var gY = svg.append('g')
                .attr('class', 'y axis')
                .call(yAxis);

            var gyLabel = gY.append('text')
                .style('fill', '#666')
                .style('font-weight', '700')
                .style('font-size', '12px')
                .attr('transform', 'translate(' + (0- margin.left / 1.25) + ',' + (height / 2) + ') rotate(-90)')
                .style('text-anchor', 'middle')
                .html(yAxisLabel || yColumn);

            // this rect acts as a layer so that zooming works anywhere in the svg. otherwise, if zoom is called on
            // just svg, zoom functionality will only work when the pointer is over a circle.
            // var zoomBox = svg.append('rect')
            //     .attr('class', 'zoom')
            //     .attr('height', innerHeight)
            //     .attr('width', innerWidth);

            var circles = svg.append('g')
                .attr('class', 'circlesContainer');

            var thresholdLines = svg.append('g').attr('class', 'thresholdLines');

            d3.select('#apply').on('click', updatePlot);


            //d3.select('#savesvg').on('click', exportSVG);

            // setup "save as svg" btn
            d3.select("#exportsvg").on("click", writeDownloadLinkV1);
            function writeDownloadLinkV1(){
                try {
                    var isFileSaverSupported = !!new Blob();
                    console.log(`is FileSaver Supported: ${isFileSaverSupported}`)
                } catch (e) {
                    alert("blob not supported");
                }
                
                
                var html = d3.select("#svg")
                    .attr("title", "test2")
                    .attr("version", 1.1)
                    .attr("xmlns", "http://www.w3.org/2000/svg")
                    .node().outerHTML; // .parentNode.innerHTML;
                
                if(!html){console.log('issue getting html: ', html);}
                var blob = new Blob([html], {type: "image/svg+xml;charset=utf-8"}); 
                //console.log('blob: ', blob);
                saveAs(blob, "exampleSVG.svg");
            }       
            // setup save as png btn
            d3.select("#exportpng").on('click', exportPNG)
            function exportPNG(){
                const output = {"name": "result.png", "width": 960, "height": 500}
                const svgElem = document.querySelector("svg")
                // const uriData = `data:image/svg+xml;base64,${btoa(svgElem.outerHTML)}` // it may fail.
                const uriData = `data:image/svg+xml;base64,${btoa(new XMLSerializer().serializeToString(svgElem))}`
                const img = new Image()
                img.src = uriData
                img.onload = () => {
                    const canvas = document.createElement("canvas");
                    [canvas.width, canvas.height] = [output.width, output.height]
                    const ctx = canvas.getContext("2d")
                    ctx.drawImage(img, 0, 0, output.width, output.height)

                    // 👇 download
                    const a = document.createElement("a")
                    const quality = 1.0 // https://developer.mozilla.org/en-US/docs/Web/API/CanvasRenderingContext2D/imageSmoothingQuality
                    a.href = canvas.toDataURL("image/png", quality)
                    a.download = output.name
                    a.append(canvas)
                    a.click()
                    a.remove()
                }
            }

            // setup reset btn
            d3.select("#resetformandplot").on('click', resetFormAndPlot)
            function resetFormAndPlot(){
                let formElement = document.getElementById('generatedplotform');

                // remove inner html
                formElement.remove();//Parent.innerHTML = '';

                // create new form and update plot
                createForm('sidebar');
                document.getElementById("defaultOpenFormTab").click();
                updatePlot();
            }

            var tooltip = d3.select("#tooltip")//("body")
                .append("div");
                //.attr('class', 'tooltip');

            var highlightedLabels = svg.append('g').attr('id', 'highlighted-labels');
                
            function tipEnter(d) {
                tooltip.style('visibility', 'visible')
                    .style('display', 'block')
                    .style('font-size', '11px')
                    .style('position', 'absolute')
                    .style('padding', '2px 7px')
                    .style('border-radius', '3px')
                    .style('background-color', '#000')
                    .style('color', '#fff')
                    .style('opacity', '0.75')
                    .style('z-index', 100)
                    .html(
                        '<p>' +
                        '<strong>' + sampleID + '</strong>: ' + d[sampleID] + '<br/>' +
                        '<strong>' + xColumn + '</strong>: ' + d3.format('.2f')(d[xColumn]) + '<br/>' +
                        '<strong>' + yColumn + '</strong>: ' + d[yColumn] +
                        '</p>'
                    )
                    .style("top", (event.pageY - 5) + "px")
                    .style("left", (event.pageX + 20) + "px");
            }

            /** 
             * @deprecated
             * moving mouse within data point won't move the tooltip (too much re-rendering)
             */
            function tipMove() {
                tooltip.style("top", (event.pageY - 5) + "px")
                    .style("left", (event.pageX + 20) + "px");
            }

            function yTickFormat(n) {
                return d3.format(".2r")(getBaseLog(10, n));
                function getBaseLog(x, y) {
                    return Math.log(y) / Math.log(x);
                }
            }

            function zoomFunction() {
                var transform = d3.zoomTransform(this);
                d3.selectAll('.dot')
                    .attr('transform', transform)
                    .attr('r', 3 / Math.sqrt(transform.k));
                gX.call(xAxis.scale(d3.event.transform.rescaleX(xScale)));
                gY.call(yAxis.scale(d3.event.transform.rescaleY(yScale)));
                svg.selectAll('.threshold')
                    .attr('transform', transform)
                    .attr('stroke-width', 1 / transform.k);
            }

            /**
             * @ deprecated
             * returned a css class name to use according to a data points position in relation to threshold lines
             *
            function circleClass(d) {
                if (d[yColumn] <= significanceThreshold && Math.abs(d[xColumn]) >= foldChangeThreshold) return 'dot sigfold';
                else if (d[yColumn] <= significanceThreshold) return 'dot sig';
                else if (Math.abs(d[xColumn]) >= foldChangeThreshold) return 'dot fold';
                else return 'dot';
            }*/

            /**
             * returns a color hex color according to a data points position in relation to threshold lines
             * hex color based on what the user has selected in 'advanced plot options' tab
             * @returns {string} hex value of color for dot in plot
             */
            function circleStyle(d) {
                // if highlighting is done using a distinct color, check if d[sampleId] is in genesToHighlightList
                // if it is, color according to highlight
                // now deprecated since highlighting is handled by parsing through data
                /*if (genesToHighlightList.includes(d[sampleID])){
                    //console.log('found data to highlight. gene name: ', d[sampleID])
                    return highlightColor;
                }*/

                // where color scheme toggle can be set
                // colors used should be able to be set by user
                //if (d[yColumn] <= significanceThreshold && Math.abs(d[xColumn]) >= foldChangeThreshold) return "#00ff00";//'green';
                if(d[yColumn] <= significanceThreshold && d[xColumn] >= foldChangeThreshold) return acceptedColorXpos;
                else if(d[yColumn] <= significanceThreshold && Math.abs(d[xColumn]) >= foldChangeThreshold) return acceptedColorXneg;
                else if (d[yColumn] <= significanceThreshold) return rejectedColorX; //"#ff0000"//'red';
                else if (Math.abs(d[xColumn]) >= foldChangeThreshold) return rejectedColorY; //'#cccccc';
                else return rejectedColorBoth; //'#000000';
            }

            /** 
             * returns a numbers for the radius of the circle based on if the data point is to be highlighted or not
             */
            function circlesSizer(d){
                // deprecated in place of parsing data in redrawPlot()
                /*if (genesToHighlightList.includes(d[sampleID])){
                    //console.log('found data to highlight. gene name: ', d[sampleID])
                    return highlightPlotDotSize;// = d3.select("#highlightcolor").property("value");
                }
                else{*/
                    return plotDotSize;
                //}
            }

            /**
             * assigns opacity to a data point on plot based on if the dot is to be highlighted or not (opcaity established in advanced tab on form)
             * @param d data point 
             * @returns float 0-1 representing level of opcacity for some dot in the plot
             */
            function circlesOpaquer(d){
                // deprecated in place of parsing data in redrawPlot()
                /*if (genesToHighlightList.includes(d[sampleID])){
                    //console.log('found data to highlight. gene name: ', d[sampleID])
                    return highlightPlotDotTransparency;// = d3.select("#highlightcolor").property("value");
                }
                else{*/
                    return plotDotTransparency;
                //}
            }

            /**
             * TODO: re-implement zoom feature, be sure to make plot redraw axes
             */
            function resetZoom() {
                var ease = d3.easePolyIn.exponent(4.0);
                svg.transition().duration(750)
                    .ease(ease)
                    .call(zoom.transform, d3.zoomIdentity);
            }

            d3.select("#export-tsv").on('click', exportTSV)
            /**
             * export tsv file with only data columns we use (include pval and padj)
             */
            function exportTSV(){
                // establish columns to include
                let columnsToInclude = new Set([dataNameCol, xColumn, yColumn, 'padj', 'pvalue']);
                /* created obj but string makes this easier
                let tsvFileObj = {}
                for(let colName of columnsToInclude){
                    tsvFileObj[colName] = [colName];
                }
                for(let dataIdx = 0; dataIdx < data.length; dataIdx++){
                    for(let colName of columnsToInclude){
                        tsvFileObj[colName].push(data[dataIdx][colName]);
                    }
                }*/
                
                // create string for tsv
                let tsvRowStrings = [];
                tsvRowStrings.push([...columnsToInclude].join('\t'));
                console.log('headers: ', tsvRowStrings)

                try{
                    // get all row data
                    for(let geneIdx=0; geneIdx < data.length; geneIdx++){
                        let curRow = [];
                        for(let colName of columnsToInclude){
                            curRow.push(data[geneIdx][colName])
                        }
                        tsvRowStrings.push(curRow.join('\t'));
                    }

                    // make blob

                    // join rows with newline char
                    let tsvFile = new File(
                        //Buffer.from(tsvRowStrings.join('\n'), 'utf-8').toString(),
                        [tsvRowStrings.join('\n')],
                        'dataFile.tsv'
                    );



                    var tsvUrl = URL.createObjectURL(tsvFile);
                    var downloadLink = document.createElement("a");
                    downloadLink.href = tsvUrl;
                    downloadLink.download = "dataFile.tsv";
                    document.body.appendChild(downloadLink);
                    //downloadLink.setAttribute('onClick', window.open(this.href,'popUpWindow'));//onclick = "window.open(this.href,'popUpWindow')"
                    downloadLink.click();
                    document.body.removeChild(downloadLink);
                }catch(err){console.log('Error creating tsv file from array of string', 'stringArr: ', tsvRowStrings, 'err: ', err)}
            }  

            /**
             * updates plot according to form values
             * steps:
             * 1. get form values
             * 2. reset relevant plot components
             * 3. redraw plot
             */
            function updatePlot() {
                // BEFORE anything else, check if aspect ratio is custom
                let isAspectRatioCustom = d3.select('#customAspectRatioToggle-selections').property('value');
                if(isAspectRatioCustom === "Custom"){
                    // set height based on static width and aspect ratio

                    let newYaspect = parseFloat(d3.select('#customAspectY').property('value'));
                    let newXaspect = parseFloat(d3.select('#customAspectX').property('value'));
                    
                    // don't allow ratios to be negative
                    if(newYaspect == 0){
                        newYaspect = 1;
                    }
                    else if (newYaspect < 0){
                        newYaspect *= -1;
                    }
                    if(newXaspect == 0){
                        newXaspect = 1;
                    }
                    else if (newXaspect < 0){
                        newXaspect *= -1;
                    }

                    
                    height = width * newYaspect / newXaspect;

                    innerWidth = width - margin.left - margin.right;
                    innerHeight = height - margin.top - margin.bottom; 


                    // reset svg data
                    //svg.remove();

                    /*svg = d3.select('#svg').append('svg')
                        .attr('id', 'svg')
                        .attr("preserveAspectRatio", "xMaxYMax slice")
                        .attr("viewBox", `0 0 ${width} ${height}`)   
                        .style('background-color', 'white')
                        .append('g')
                        .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

                    svg.append('defs').append('clipPath')
                        .attr('id', 'clip')
                        .append('rect')
                        .attr('height', innerHeight)
                        .attr('width', innerWidth);
                    */

                    // try editing attr instead of re making svg
                    d3.select('#svg').attr("viewBox", `0 0 ${width} ${height}`);
                    d3.select('#clip').attr('height', innerHeight).attr('width', innerWidth);
                }
                else{
                    height = 500;
                    innerWidth = width - margin.left - margin.right;
                    innerHeight = height - margin.top - margin.bottom; 


                    // try editing attr instead of re making svg
                    d3.select('#svg').attr("viewBox", `0 0 ${width} ${height}`);
                    d3.select('#clip').attr('height', innerHeight).attr('width', innerWidth);
                }


                /* 
                    var ease = d3.easePolyIn.exponent(4.0);
                    svg.transition().duration(750)
                        .ease(ease)
                        .call(zoom.transform, d3.zoomIdentity);
                

                console.log(`width: ${width}, height: ${height}, innerW: ${innerWidth}, innerH: ${innerHeight}`)
                */

                // ************* get values from form *************

                plotDotSize = parseInt(d3.select("#plotdotsize").property("value"));
                plotDotTransparency = parseFloat(d3.select("#plotdotopacity").property("value"));
                foldChangeThreshold = parseFloat(d3.select("#foldchange").property("value"));
                significanceThreshold = parseFloat(d3.select("#fdrthreshold").property("value"));
                plotTitle = d3.select('#plottitle').property("value");
                xAxisLabel = d3.select("#plotxaxis").property("value");
                yAxisLabel = d3.select("#plotyaxis").property("value");
                acceptedColorXpos = d3.select("#acceptedxcolorpos").property("value");
                acceptedColorXneg = d3.select("#acceptedxcolorneg").property("value");
                rejectedColorX = d3.select("#rejectedxcolor").property("value");
                rejectedColorY = d3.select("#rejectedycolor").property("value");
                rejectedColorBoth = d3.select("#rejectedbothcolor").property("value");
                highlightColor = d3.select("#highlightcolor").property("value");
                highlightPlotDotSize = parseInt(d3.select("#highlightedplotdotsize").property("value"));
                highlightPlotDotTransparency = parseFloat(d3.select("#highlightedplotdotopacity").property("value"));
                genesToHighlightList = [];
                highlightTopXBy = d3.select('#highlighttopn-selections-by').property('value'); // by for what value to use when sorting
                highlightTopXFilter = d3.select('#highlighttopn-selections-filter').property('value');
                highlightTopXGenes = parseFloat(d3.select('#highlighttopn').property("value"));
                d3.select('#highlightlistford3').each(function (p, j) {
                    //console.log(`genes to highlight. p: ${p}, j: ${j}`);
                    d3.select(this)
                        .selectAll("li")
                        .each(function (ele, ctx) {
                            let tempD = d3.select(this);
                            //console.log('tempD.node().value', tempD.node().value, tempD.property('innerText'), tempD);
                            genesToHighlightList.push(tempD.property('innerText'));
                        })
                });
                //console.log(`genes to highlight: ${genesToHighlightList}`);
                let customAxisX = d3.select('#customxrangetoggle-selections').property('value');
                let customAxisY = d3.select('#customyrangetoggle-selections').property('value');

                // set highlightTopxBY based on value (will include By: )
                highlightTopXBy = highlightTopXBy.split(' ')[1];
                //highlightTopXFilter = highlightTopXFilter.split(' ')[1];
                //console.log('highlightTopXBy: ', highlightTopXBy);
                console.log('highlightTopXFilter: ', highlightTopXFilter);
                // reset the plot
                highlightedLabels.remove()
                thresholdLines.remove();
                circles.remove();
                curPlotTitle.remove();
                gxLabel.remove();
                gyLabel.remove();
                gY.remove();
                gX.remove();


                // *** handle when axis are custom ***

                if(customAxisX === 'Custom'){
                    console.log("make custom x axis range");
                    // set x scale
                    let newXmin = parseFloat(d3.select('#customxrangemin').property('value'));
                    let newXmax = parseFloat(d3.select('#customxrangemax').property('value'));
                    xScale.range([0, innerWidth])
                        .domain([newXmin, newXmax])//d3.extent(data, function (d) { return d[xColumn]; }))
                        .nice();
                }
                else{
                    xScale.range([0, innerWidth])
                        .domain(d3.extent(data, function (d) { return d[xColumn]; }))
                        .nice();
                    
                }
                if(customAxisY === 'Custom'){
                    console.log("make custom y axis range");
                    let newYmin = parseFloat(`1e${d3.select('#customyrangemin').property('value')}`);
                    let newYmax = 1;//parseFloat(d3.select('#customyrangemax').property('value'));
                    yScale.range([0, innerHeight])
                        .domain([newYmin, newYmax])//d3.extent(data, function (d) { return d[yColumn]; }))
                        .nice();
                }
                else{
                    let yScaleDomain;
                    try {
                        yScaleDomain = d3.extent(data, function (d) { return d[yColumn]; })
                        console.log('default yscale domain: ', yScaleDomain);
                        // if lower bound is 0, set to usable number
                        if(yScaleDomain[0] == 0) {
                            yScaleDomain[0] = parseFloat('1e-308');
                        }
                    } catch (err) {
                        console.log('issue creating yscale default domain. err: ', err);
                        console.log('yscale domain BEFORE altered set: ', yScaleDomain);
                        yScaleDomain = [parseFloat(`1e-308`), 1];
                        console.log('yscale domain AFTER altered set: ', yScaleDomain);
                    }
                    yScale.range([0, innerHeight])
                        .domain(yScaleDomain)
                        .nice();
                }

                console.log('xScale domain: ', xScale.domain())
                console.log('yScale domain: ', yScale.domain())
                console.log(typeof yScale.domain()[0], yScale.domain()[0])

                // ************* redraw the plot *************
                xAxis = d3.axisBottom(xScale);
                yAxis = d3.axisLeft(yScale)
                    .ticks(5)
                    .tickFormat(yTickFormat);

                gX = svg.append('g')
                    .attr('class', 'x axis')
                    .attr('transform', 'translate(0,' + innerHeight + ')')
                    .call(xAxis);
                gY = svg.append('g')
                    .attr('class', 'y axis')
                    .call(yAxis);

                curPlotTitle = svg.append("text")
                    .attr('transform', 'translate(' + width / 2 + ',' + (margin.top - 25) + ')')
                    .attr("text-anchor", "middle")  
                    .style("font-size", "16px") 
                    .style("text-decoration", "underline")  
                    .text(plotTitle);

                let xLabelText = `Log2 ${xAxisLabel || xColumn}`;
                gxLabel = gX.append('text')
                    .style('fill', '#666')
                    .style('font-weight', '700')
                    .style('font-size', '12px')
                    .attr('transform', 'translate(' + width / 2 + ',' + (margin.bottom - 6) + ')')
                    .attr('text-anchor', 'middle')
                    .html(xLabelText);
                
                let yLabelText = `Log10 ${yAxisLabel || yColumn}`;
                gyLabel = gY.append('text')
                    .style('fill', '#666')
                    .style('font-weight', '700')
                    .style('font-size', '12px')
                    .attr('transform', 'translate(' + (0 - margin.left / 1.25) + ',' + (height / 2) + ') rotate(-90)')
                    .style('text-anchor', 'middle')
                    .html(yLabelText);

                circles = svg.append('g')
                    .attr('class', 'circlesContainer');

                circles.selectAll(".dot")
                    .data(data)
                    .enter().append('circle')
                    .attr('id', function(d) {return `plotdot${d[sampleID]}`} )
                    .attr('r', circlesSizer) 
                    .attr('cx', function (d) { return xScale(d[xColumn]); })
                    .attr('cy', function (d) { return yScale(d[yColumn]); })
                    .style('fill', circleStyle)
                    .style('opacity', circlesOpaquer)
                    .on('mouseenter', tipEnter)
                    //.on("mousemove", tipMove)
                    .on('mouseleave', function (d) {
                        return tooltip.style('visibility', 'hidden');
                    });

                thresholdLines = svg.append('g').attr('class', 'thresholdLines');

                // add horizontal line at significance threshold
                thresholdLines.append("svg:line")
                    .style('stroke', '#333')
                    .style('stroke-dasharray', '5px 10px')
                    .style('stroke-opacity', '0.35')
                    .attr("x1", 0)
                    .attr("x2", innerWidth)
                    .attr("y1", yScale(significanceThreshold))
                    .attr("y2", yScale(significanceThreshold));

                // add vertical line(s) at fold-change threshold (and negative fold-change)
                [foldChangeThreshold, -1 * foldChangeThreshold].forEach(function (threshold) {
                    thresholdLines.append("svg:line")
                        .style('stroke', '#333')
                        .style('stroke-dasharray', '5px 10px')
                        .style('stroke-opacity', '0.35')
                        .attr("x1", xScale(threshold))
                        .attr("x2", xScale(threshold))
                        .attr("y1", 0)
                        .attr("y2", innerHeight);
                });

                /** list of objects in form { geneId: <val>, pvalue: <val>, padj: <val>,  negLOG10pval: <val>, negLOG10padj: <val> } ] */
                let pvalueMap = [] 
                highlightedLabels = svg.append('g').attr('id', 'highlighted-labels');
                // handle named (given) highlighted plot points (and build temp array for sorting and filtering data)
                for (let geneI=0; geneI<data.length; geneI++) {
                    let dataName = data[geneI][sampleID];
                    let dataValX = data[geneI][xColumn]; 
                    let dataValY = data[geneI][yColumn];
                    if (geneListToHighlight.includes(dataName)){
                        console.log(`n: ${dataName}, x: ${dataValX}, y: ${dataValY}`);                        
                        
                        // remove and redraw circles (remove to put highlighted points on top of others)
                        d3.select(`#plotdot${dataName}`).remove();

                        circles.append('circle')
                            //.data({sampleID: dataName, yColumn: dataValY, xColumn: dataValX})
                            //.enter()
                            .attr('id', function(d) {return `plotdot${dataName}`} )
                            .attr('r', highlightPlotDotSize) 
                            .attr('cx', function (d) { return xScale(dataValX); })
                            .attr('cy', function (d) { return yScale(dataValY); })
                            .style('fill', highlightColor)
                            .style('opacity', highlightPlotDotTransparency)
                            .style('z-index', 100)
                            .on('mouseenter', function(){return tipEnter({sampleID: dataName, yColumn: dataValY, xColumn: dataValX});})// function(d){ return tipEnter(d)})
                            //.on("mousemove", tipMove)
                            .on('mouseleave', function (d) {
                                return tooltip.style('visibility', 'hidden');
                            });
                        // tooltip not able to be applied

                        
                        
                        

                        // draw tooltip (needs to be apart of the svg if its to be included in the exported image)
                        highlightedLabels.append('g')
                            .attr('id', `highlighttooltip${dataName}`)
                            .attr('transform', `translate(${xScale(dataValX) - 5},${yScale(dataValY) + 20})`)
                            //.style('position', 'absolute')
                            //.style('padding', '2px 7px')
                            //.style('border-radius', '3px')
                            //.style('fill', '#000')//.style('background-color', '#000')
                            //.style('opacity', '0.75')
                            .style('z-index', 100)
                            .append('text')
                            .text(`${dataName}`) //`${sampleID}: ${dataName}`
                            .style('font-size', '11px')
                            //.style('color', '#fff')
                            //*/
                            //.attr("text-anchor", "middle")  
                            //.style("text-decoration", "underline")
                        
                        console.log('test')
                    }
                    
                    // TODO: save off geneID and adj. pval for sorting and getting top X to highlight
                    let tempPvalObj = {
                        name: dataName,
                        pvalue: data[geneI]['pvalue'],
                        padj: data[geneI]['padj'],
                        //negLOG10pval: parseFloat(data[geneI]["'-LOG10(pval)'"]),
                        //negLOG10padj: data[geneI]["'-LOG10(padj)'"],

                        // save x and y for redrawing circle
                        xVal: dataValX,
                        yVal: dataValY
                    }
                    pvalueMap.push(tempPvalObj);
                }
                // testing how to access pvals and sort to find genes to highlight
                // console.log(`data cols: ${data['columns']}`);
                // console.log('pvalue map: ', pvalueMap);
                // console.log(data[0]);

                // before sorting, filter pvalueMap by threshold
                pvalueMap = pvalueMap.filter((dataPoint) => {
                    let thresholdStatus;
                    if(dataPoint.xVal >= foldChangeThreshold){thresholdStatus = "pos";}
                    else if(Math.abs(dataPoint.xVal) >= foldChangeThreshold){thresholdStatus = 'neg'}
                    else{thresholdStatus = "fail"}

                    if(highlightTopXFilter.includes('all')){
                        return true;
                    }
                    else if(highlightTopXFilter.includes('both') && (thresholdStatus == "pos" || thresholdStatus == "neg")){
                        return true;
                    }
                    else if(highlightTopXFilter.includes('neither') && thresholdStatus == "fail"){
                        return true;
                    }
                    else if(highlightTopXFilter.includes('negative') && thresholdStatus == "neg"){
                        return true;
                    }
                    else if(highlightTopXFilter.includes('positive') && thresholdStatus == "pos"){
                        return true;
                    }
                    return false
                })

                // sort pvalueMap by selected metric (pvalue vs padj as options)
                pvalueMap.sort((a, b) => (a[highlightTopXBy] - b[highlightTopXBy]));//((a, b) => (a[highlightTopXBy] > b[highlightTopXBy])) ? 1 : 0); //((a.pvalue > b.pvalue) ? 1 : 0); //(ttt, (a, b) => (a[highlightTopXBy] - b[highlightTopXBy]));
                //console.log(pvalueMap.slice(0, 10))
                if(!!highlightTopXGenes && highlightTopXGenes > 0){
                    for(let i=0;i<highlightTopXGenes;i++){ //for top X genes (by pval)
                        if(!geneListToHighlight.includes(pvalueMap[i].name)){ // if gene not already manaully highlighted
                            let dataName = pvalueMap[i].name,
                                dataValX = pvalueMap[i].xVal,
                                dataValY = pvalueMap[i].yVal;

                            // remove and redraw circles (remove to put highlighted points on top of others)
                            d3.select(`#plotdot${dataName}`).remove();

                            circles.append('circle')
                                .attr('id', function(d) {return `plotdot${dataName}`} )
                                .attr('r', highlightPlotDotSize) 
                                .attr('cx', function (d) { return xScale(dataValX); })
                                .attr('cy', function (d) { return yScale(dataValY); })
                                .style('fill', highlightColor)
                                .style('opacity', highlightPlotDotTransparency)
                                .style('z-index', 100)
                                .on('mouseenter', function(){return tipEnter({sampleID: dataName, yColumn: dataValY, xColumn: dataValX});})// function(d){ return tipEnter(d)})
                                .on('mouseleave', function (d) {
                                    return tooltip.style('visibility', 'hidden');
                                });
                            // tooltip not able to be applied
                            

                            // draw tooltip (needs to be apart of the svg if its to be included in the exported image)
                            highlightedLabels.append('g')
                                .attr('id', `highlighttooltip${dataName}`)
                                .attr('transform', `translate(${xScale(dataValX) - 5},${yScale(dataValY) + 20})`)
                                .style('z-index', 100)
                                .append('text')
                                .text(` ${dataName}`)
                                .style('font-size', '11px')
                        
                        }
                    }
                }

            } 


            
            /**
             * add funtion to d3 to support adding labels for highlighted genes
             *
            d3.force_labels = function force_labels() {    
                var labels = d3.layout.force();
                  
                // Update the position of the anchor based on the center of bounding box
                function updateAnchor() {
                  if (!labels.selection) return;
                  labels.selection.each(function(d) {
                    var bbox = this.getBBox(),
                        x = bbox.x + bbox.width / 2,
                        y = bbox.y + bbox.height / 2;
            
                    d.anchorPos.x = x;
                    d.anchorPos.y = y;
                   
                    // If a label position does not exist, set it to be the anchor position 
                    if (d.labelPos.x == null) {
                      d.labelPos.x = x;
                      d.labelPos.y = y;
                    }
                  });
                }
                
                //The anchor position should be updated on each tick
                labels.on("tick.labels", updateAnchor);
                
                // This updates all nodes/links - retaining any previous labelPos on updated nodes
                labels.update = function(sel) {
                  labels.selection = sel;
                  var nodes = [], links = [];
                  sel[0].forEach(function(d) {    
                    if(d && d.__data__) {
                      var tempData = d.__data__;
                      
                      if (!d.labelPos) d.labelPos = {fixed: false};
                      if (!d.anchorPos) d.anchorPos = {fixed: true};
                      
                      // Place position objects in __data__ to make them available through 
                      // d.labelPos/d.anchorPos for different elements
                      tempData.labelPos = d.labelPos;
                      tempData.anchorPos = d.anchorPos;
                      
                      links.push({target: d.anchorPos, source: d.labelPos});
                      nodes.push(d.anchorPos);
                      nodes.push(d.labelPos);
                    }
                  });
                  labels
                      .stop()
                      .nodes(nodes)
                      .links(links);
                  updateAnchor();
                  labels.start();
                };
                return labels;
            }*/
            

            updatePlot();
        });
    }

    // all functions below are used to set aspects of the plot
    /*chart.genesList = function (value) {
        if (!arguments.length) return genesList;
        genesList = value;
        return chart;
    };*/

    chart.width = function (value) {
        if (!arguments.length) return width;
        width = value;
        return chart;
    };

    chart.height = function (value) {
        if (!arguments.length) return height;
        height = value;
        return chart;
    };

    chart.margin = function (value) {
        if (!arguments.length) return margin;
        margin = value;
        return chart;
    };

    chart.xColumn = function (value) {
        if (!arguments.length) return xColumn;
        xColumn = value;
        return chart;
    };

    chart.yColumn = function (value) {
        if (!arguments.length) return yColumn;
        yColumn = value;
        return chart;
    };

    chart.xAxisLabel = function (value) {
        if (!arguments.length) return xAxisLabel;
        xAxisLabel = value;
        return chart;
    };

    chart.yAxisLabel = function (value) {
        if (!arguments.length) return yAxisLabel;
        yAxisLabel = value;
        return chart;
    };

    chart.xAxisLabelOffset = function (value) {
        if (!arguments.length) return xAxisLabelOffset;
        xAxisLabelOffset = value;
        return chart;
    };

    chart.yAxisLabelOffset = function (value) {
        if (!arguments.length) return yAxisLabelOffset;
        yAxisLabelOffset = value;
        return chart;
    };

    chart.xTicks = function (value) {
        if (!arguments.length) return xTicks;
        xTicks = value;
        return chart;
    };

    chart.yTicks = function (value) {
        if (!arguments.length) return yTicks;
        yTicks = value;
        return chart;
    };

    chart.significanceThreshold = function (value) {
        if (!arguments.length) return significanceThreshold;
        significanceThreshold = value;
        return chart;
    };

    chart.foldChangeThreshold = function (value) {
        if (!arguments.length) return foldChangeThreshold;
        foldChangeThreshold = value;
        return chart;
    };

    chart.colorRange = function (value) {
        if (!arguments.length) return colorRange;
        colorRange = value;
        return chart;
    };

    chart.sampleID = function (value) {
        if (!arguments.length) return sampleID;
        sampleID = value;
        return chart;
    };
    
    chart.plotTitle = function (value) {
        if (!arguments.length) return plotTitle;
        plotTitle = value;
        return chart;
    };

    return chart;
}

/**
 * create's the interactive pages form
 * 
 * attaches bootstrap card to form 
 * attaches tabs to card header (default tab given id so it's open by deafult)
 * attaches tab content to divs within card-content with ID's correlating to their tab
 */
function createForm(formDivId){
    /** object for structuring form content
     * 
     * object follows following structure
     * 
     * formObj: {
            tabNames: array of names of tabs in form. also name of fields in object structured for that tabs form,
            <tab_name_1>: {
                displayName: string to display for tabs title (clickable for swapping between tabs in form),
                htmlId: id of the div to contain this tabs forms,
                formControls: {
                    order: list of different controls needed in form (in order they are to exist in on form) (also name of fields in this object),
                    <order_name_1>: {
                        label: string to be displayed on form for differentiating inputs,
                        inputType: type to be injeted into html as input type,
                        value: value for the form to initially hold,
                        htmlId: id of the input form,
                        (optional) extraAddons: {
                            (optional) divs: {
                                order: [] list of names following form control structure for divs to be added to this formControls input-group,
                                <div_name_1>: {
                                    (optional) label: label to be give to div

                                },
                                ...
                            },
                            (optional) styles: {} each field is name of style, value is value for that style,
                            (optional) options: [{ //list of selectable input options
                                (optional) label: label prepended to this dropdown selections input group,
                                (optional) selectorID: identifier for when this list has more than 1 selectable dropdown,
                                type: determines selectability ('single-select', 'multi-select')
                                options: [] list of options to display

                            }, ...],
                            (optional) toAppend: {} // formatted the same as form controls but for adding other html elements to a form control
                        }
                    },
                    ...
                },
            },
            ...
        }

        NOTES: 
         - formControls with inputTypes of "none" are meant to be created as divs (within a separate input group on that tabs form)
     */
    var webFormContent = {
        tabNames: ["main", "advanced", "color", "highlight"],
        main: {
            displayName: "Basic Plot Options",
            htmlId: 'mainform',
            formControls: {
                order: ['title', 'ylabel', 'xlabel', 'fdrthreshold', 'foldchangethreshold'],
                title: {
                    label: 'Plot Title: ',
                    inputType: 'text',
                    value: 'Plot Title',
                    htmlId: 'plottitle',
                },
                ylabel: {
                    label: 'Y Axis Label: ',
                    inputType: 'text',
                    value: 'False Discovery Rate',//'-log<tspan baseline-shift="sub">10</tspan>False Discovery Rate',
                    htmlId: 'plotyaxis',
                },
                xlabel: {
                    label: 'X Axis Label: ',
                    inputType: 'text',
                    value: 'Fold-change',//'log<tspan baseline-shift="sub">2</tspan>Fold-change',
                    htmlId: 'plotxaxis',
                },
                fdrthreshold: {
                    label: '-Log10 FDR threshold: ',
                    inputType: 'text',
                    value: '0.05',
                    htmlId: 'fdrthreshold',
                },
                foldchangethreshold: {
                    label: 'log2 Fold-change threshold: ',
                    inputType: 'text',
                    value: '1',
                    htmlId: 'foldchange',
                }
            }
        },
        advanced: {
            displayName: "Advanced Plot Options",
            htmlId: "advancedoptions",
            formControls: {
                /** Requirements
                 *  toggle to tell updatePlot method to use custom x and/or y ranges (separate toggle for each?)
                 *  4 input boxes for min and max (for x&y)
                 * 
                 * changes to volcanoPlot needed:
                 *  updatePlot should check if x or y toggle is checked (before circles anything is redrawn)
                 *  if checked, get new ranges and set x/y Scales .domain() based on this range
                 */
                order: ['customXtoggle', 'customXminRange', 'customXmaxRange', 'customYtoggle', 'customYminRange', /*'customYmaxRange',*/ 'customAspectRatioToggle', "customAspectX", "customAspectY"],
                customXtoggle: {
                    label: "Use custom range for x-axis", 
                    inputType: 'none',
                    htmlId: 'customxrangetoggle',
                    extraAddons: {
                        options: [
                            {
                                type: "single-select", // only 1 of the options may be selected at once
                                options: ["Default", "Custom"]
                            }
                        ]
                    }
                },
                customXminRange: {
                    label: "x-axis (log2) min", 
                    inputType: 'number',
                    value: "-8",
                    htmlId: 'customxrangemin',
                },
                customXmaxRange: {
                    label: "x-axis (log2) max", 
                    inputType: 'number',
                    value: "10",
                    htmlId: 'customxrangemax',
                },
                customYtoggle: {
                    label: "Use custom range for y-axis", 
                    inputType: 'none',
                    htmlId: 'customyrangetoggle',
                    extraAddons: {
                        options: [
                            {
                                type: "single-select", // only 1 of the options may be selected at once
                                options: ["Default", "Custom"]
                            }
                        ]
                    }
                },
                customYminRange: {
                    label: "y-axis (-log10) min", 
                    inputType: 'number',
                    value: "-308",
                    htmlId: 'customyrangemin',
                },
                // customYmaxRange: {
                //     label: "y-axis (-log10) max (should always be 1)", 
                //     inputType: 'number',
                //     value: "1",
                //     htmlId: 'customyrangemax',
                // },
                
                customAspectRatioToggle: {
                    label: "Use custom aspect ratio for plot (width-INT : height-INT)", 
                    inputType: 'none',
                    htmlId: 'customAspectRatioToggle',
                    extraAddons: {
                        options: [
                            {
                                type: "single-select", // only 1 of the options may be selected at once
                                options: ["Default", "Custom"]
                            }
                        ]
                    }
                },
                customAspectX: {
                    label: "Custom X Aspect",
                    inputType: "number",
                    value: "6",
                    htmlId: "customAspectX"
                },
                customAspectY: {
                    label: "Custom Y Aspect",
                    inputType: "number",
                    value: "4",
                    htmlId: "customAspectY"
                }
            }
        },
        color: {
            displayName: "Color Options",
            htmlId: 'colorform',
            formControls: {
                order: ['plotdotopacity', 'plotdotsize', 'acceptedxcolorpos', 'acceptedxcolorneg', 'rejectedxcolor', 'rejectedycolor', 'rejectedbothcolor'],  
                plotdotopacity: {
                    label: 'data point transparency (less than 1): ',
                    inputType: 'text',
                    value: '0.4',
                    htmlId: 'plotdotopacity',
                },
                plotdotsize: {
                    label: 'data point size (1-100): ',
                    inputType: 'number',
                    value: '3',
                    htmlId: 'plotdotsize',
                    extraAddons: {
                        styles: {
                            min: "1",
                            max: "100"
                        }
                    }
                },
                acceptedxcolorpos: {
                    label: 'Accepted (positive X) color: ',
                    inputType: 'color',
                    value: '#00ff00',
                    htmlId: 'acceptedxcolorpos'
                },
                acceptedxcolorneg: {
                    label: 'Accepted (negative X) color: ',
                    inputType: 'color',
                    value: '#00ff00',
                    htmlId: 'acceptedxcolorneg'
                },
                rejectedxcolor: {
                    label: 'Rejected (X only) color: ',
                    inputType: 'color',
                    value: '#ff0000',
                    htmlId: 'rejectedxcolor'
                },
                rejectedycolor: {
                    label: 'Rejected (Y only) color: ',
                    inputType: 'color',
                    value: '#cccccc',
                    htmlId: 'rejectedycolor'
                },
                rejectedbothcolor: {
                    label: 'Rejected (X and Y) color: ',
                    inputType: 'color',
                    value: '#000000',
                    htmlId: 'rejectedbothcolor'
                },
            }
        },
        highlight: {
            displayName: "Highlight Options",
            htmlId: 'highlightform',
            formControls: {
                order: ['highlightedplotdotopacity', 'highlightedplotdotsize', 'highlightcolor', 'highlighttopn', 'genehighlight', 'selectedGeneList'],
                highlightedplotdotopacity: {
                    label: 'highlighted data point transparency (less than 1): ',
                    inputType: 'text',
                    value: '0.4',
                    htmlId: 'highlightedplotdotopacity',
                },
                highlightedplotdotsize: {
                    label: 'highlighted data point size (1-100): ',
                    inputType: 'number',
                    value: '3',
                    htmlId: 'highlightedplotdotsize',
                    extraAddons: {
                        styles: {
                            min: "1",
                            max: "100"
                        }
                    }
                },
                highlightcolor: {
                    label: 'Highlighted data point color: ',
                    inputType: 'color',
                    value: '#abcdef',
                    htmlId: 'highlightcolor',
                },
                highlighttopn: {
                    label: 'Highlight top X genes: ',
                    inputType: 'number',
                    value: '0',
                    htmlId: 'highlighttopn',
                    extraAddons: {
                        styles: {
                            min: "1",
                            max: "100"
                        },
                        options: [
                            {
                                label: "by: ",
                                selectorID: "by",
                                type: "single-select", 
                                options: ["pvalue", "padj"]
                            },
                            {
                                label: "include: ",
                                selectorID: "filter",
                                type: "single-select", 
                                options: ["all", "positive of X threshold only", "negative of X threshold only", "both (pos/neg) of X threshold", "neither (pos/neg) of X threshold"]
                            }
                        ]
                    }
                },
                genehighlight: {
                    label: 'Enter gene name(s) (click name on right): ',
                    inputType: 'text',
                    value: '',
                    htmlId: 'q',
                    extraAddons: {
                        styles: {
                            name: "q",
                            onKeyUp: "showResults(this.value)"
                        },
                        divs: {
                            order: ['searchResult'],
                            searchResult: {
                                htmlId: 'genesearchresult'
                            }
                        }
                    }
                },
                selectedGeneList: {
                    label: "Selected Genes to highlight: ",
                    htmlId: "genehighlightlist",
                    inputType: "none",
                    value: ""
                }
            }
        },

    }

    var plotForm = d3.select(`#${formDivId}`).append('div') //'#sidebar').append('div')
        .attr('class', 'card')
        .attr('id', 'generatedplotform');
    

    /**
     * adds tabs for form sections
     * default tab should be named "main" in the {webFormContent} object, so it can be given an ID that will auto select it when the page initially loads 
     */
    function addFormTabs(){
        var tabHeader = plotForm.append('div')
            .attr('class', 'card-header') //card-header // row tab-header
            .append('div')
            .attr('class', 'tab row');
        
        // add buttons for each tab
        for( let tabName of webFormContent.tabNames){
            // main tab gets an id so that it can be open by default (also formats form to only show 1 tab at a time)
            if(tabName == 'main'){
                tabHeader.append('div')
                    .attr('class', 'tablinks btn btn-link')
                    .attr('onclick', `openFormTab(event, '${webFormContent[tabName].htmlId}')`)
                    .attr('id', 'defaultOpenFormTab')
                    .text(webFormContent[tabName].displayName);
            }
            else{
                tabHeader.append('div')
                    .attr('class', 'tablinks btn btn-link')
                    .attr('onclick', `openFormTab(event, '${webFormContent[tabName].htmlId}')`)
                    .text(webFormContent[tabName].displayName);
            }
        }
    }

    /**
     * includes add on divs/styles to some input control
     * @param inputGroup d3 selection of the group containing some input control (for adding divs)
     * @param inputControl d3 selection of the form-control that is getting extra styling
     * @param addOns object of styles/divs to include for some form-control
     */
    function includeAddons(inputGroup, inputControl, addOns, formControlHtmlID){
        if('divs' in addOns){
            // add divs to tabForm
            var divs = addOns.divs.order;
            for(let divToAdd of divs){
                inputGroup.append('div').attr('id', addOns.divs[divToAdd].htmlId).attr('class', 'input-group-append');
            }
        }
        if('styles' in addOns){
            //add attributes to inputControl
            for(let addOnKey of Object.keys(addOns['styles'])){
                inputControl.attr(addOnKey, addOns['styles'][addOnKey]);
            }
            
        }
        if('options' in addOns){
            // for each selectable option to add
            for(let selectableOptionObj of addOns['options']){
                // add selectable options based on type
                //let selectableOptionObj = addOns['options'];
                let promptDiv = d3.select(`#${formControlHtmlID}-prompt`);
                
                /*switch(selectableOptionObj['type']){
                    case "single-select":*/
                if(selectableOptionObj['type'] == "single-select"){
                        // for single select, 1 button with dropdown options created. text for button set based on active selection
                        let selections;

                        // before adding selectable options, check if selections are labeled and add to value
                        if(Object.keys(selectableOptionObj).includes('label')){
                            let newOptions = []
                            for(let selectableOptionText of selectableOptionObj['options']){
                                newOptions.push(`${selectableOptionObj['label']}${selectableOptionText}`);
                            }
                            selectableOptionObj['options'] = newOptions;
                        }

                        // if selectorID given, include in htmlID
                        if(Object.keys(selectableOptionObj).includes('selectorID')){
                            let distinctHtmlId = selectableOptionObj['selectorID'];
                            selections = promptDiv.append('select')
                                .attr('class', 'select')
                                .attr('id', `${formControlHtmlID}-selections-${distinctHtmlId}`)
                            // append label to option value
                            let selectionOptions = selections.selectAll('option')
                                .data(selectableOptionObj.options).enter()
                                .append('option')
                                .attr('id', `${formControlHtmlID}-option-${distinctHtmlId}`)
                                .text(function (d) { return d; })
                        }
                        else{
                            selections = promptDiv.append('select')
                                .attr('class', 'select')
                                .attr('id', `${formControlHtmlID}-selections`)
                                /*.on('change', function(){
                                    let selectVal = d3.select(this).property('value');
                                    console.log('selected val: ', selectVal)
                                    d3.select(this).text(selectVal);
                                });*/
                        
                            let selectionOptions = selections.selectAll('option')
                                .data(selectableOptionObj.options).enter()
                                .append('option')
                                .attr('id', `${formControlHtmlID}-option`)
                                .text(function (d) { return d; })
                                /*.on('click', function(){
                                    let htmlIdForSelector = d3.select(this).property('id').split("-")[0]
                                    let curValue = d3.select(this).property('value');
                                    console.log(`htmlID: ${htmlIdForSelector}`, `value clicked: ${curValue}`);
                                    let parentSelector = d3.select(`#${htmlIdForSelector}-selections`)
                                    parentSelector.text(curValue)
                                });*/
                        }
                        //break;
                    /*case "multi-select":

                        break;
                    case "toggle": 

                        break;
                    default:
                        break;*/
                    // elif multi-select
                    // elif ...
                }

            }
        }
    }

    // adds form controls for a given tab
    function addFormControlsForTabs(tabName){
        // add div for the content of a tab
        let tabForm = plotForm.append('div')
            .attr('id', webFormContent[tabName].htmlId)
            .attr('class', 'tabcontent card-body')
            .append('form');

        // get formControls for this tab
        let control = webFormContent[tabName].formControls;
        // add each form control to the tab
        for(let formControlName of control.order){
            // add label for form control
            tabForm.append('div')
                .attr('id', `${control[formControlName].htmlId}-prompt`)
                .append('p')
                //.attr('class', 'input-group-prepend')
                //.append('span')
                //.attr('class', 'input-group-text')
                .text(control[formControlName].label);
            var inputGroup = tabForm.append('div')
            // if form control is for autocompleting sample searching, include dropdown
            if(control[formControlName].inputType =="text" && control[formControlName].htmlId == "q"){ 
                inputGroup.attr('class', 'input-group mb-3 dropdown');
            }
            else{
                inputGroup.attr('class', 'input-group mb-3');
            }

            if(control[formControlName].inputType !== "none"){
                /*
                // add label for form control
                tabForm.append('div')
                    .attr('class', 'input-group-prepend')
                    .append('span')
                    .attr('class', 'input-group-text')
                    .text(control[formControlName].label);
                */
                // add input
                var inputControl = inputGroup.append('input')
                    .attr('class', 'form-control')
                    .attr('type', control[formControlName].inputType)
                    .attr('id', control[formControlName].htmlId)
                    .attr('value', control[formControlName].value)

                
            }
            else{
                // in this case, the form-control we are adding isn't for inputs, but for displaying and holding data based on inputs so that it's gettable when redrawing the plot
                
                /*
                // add label for form control
                tabForm.append('div')
                    .attr('class', 'input-group-prepend')
                    .append('span')
                    .attr('class', 'input-group-text')
                    .text(control[formControlName].label);
                */

                // add div itself with just ID
                inputGroup.append('div').attr('id', control[formControlName].htmlId);
            }
            // if there are addOns, include them
            if('extraAddons' in control[formControlName] ){
                includeAddons(inputGroup, inputControl, control[formControlName]['extraAddons'], control[formControlName].htmlId);
            }
        }
    }


    addFormTabs();

    for(let tabName of webFormContent.tabNames){
        addFormControlsForTabs(tabName);
    }
}


// for listing genes and making them searchable, these should be global and accesible by both the plot and the form
var geneList = [];
var geneListToHighlight = [];

/**
 * ? should the column names come as a list of strings for selectable options ?
 * @param {string} fileUrl url or path of file (path if used as "linkList" VP for cwl output, url if accessing data file via satellite file-manager)
 * @param {string} plotDivID id of div in html that will contain the form for editing the plot
 * @param {string} formDivID id of the div in html that will contain the plot
 * @param {string} xColName: name of column in tsv file that correlates to the plots x axis
 * @param {string} yColName: name of column in tsv file that correlates to the plots y axis
 * @param {string} dataNameCol: name of column in tsv file that correlates to the data point names (gene id’s)
 */
function initVisualPlugin(f, plotDivID, formDivID, xColName, yColName, dataNameCol){
    // create form
    createForm(formDivID)// formDivID

    // make sure main form content is open (if not here then all form tabs populate the form)
    document.getElementById("defaultOpenFormTab").click();
    
    // create volcano plot

    var volPlot = volcanoPlot()
        //.foldChangeThreshold(1.0)
        .sampleID(dataNameCol)
        .xColumn(xColName) 
        .yColumn(yColName);

    
    d3.tsv(f, parser, function (error, data) {
        if (error) console.log(error);

        // generate list of genes for selecting highlights
        for (let geneI=0; geneI<data.length; geneI++) {
            geneList.push(data[geneI][dataNameCol]);
        }
        // console.log('genelist len: ', geneList.length);
        d3.select(`#${plotDivID}`)//'#chart') //plotDivID
            .data([data])
            .call(volPlot);
    });



}

// ***** helper functions used by form and plot ********

// row parser to convert key values into numbers if possible
function parser(d) {
    for (var key in d) {
        if (d.hasOwnProperty(key)) {
            d[key] = numberParser(d[key]);
        }
    }
    return d;
}

// function to turn string into number if possible
function numberParser(value) {
    return (+value) ? +value : value;
}


// other helpers
// opens clicked tab of the form
function openFormTab(evt, tabName) {
    // Declare all variables
    var i, tabcontent, tablinks;

    // Get all elements with class="tabcontent" and hide them
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
    }

    // Get all elements with class="tablinks" and remove the class "active"
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" tablinks-active", "");
    }

    // Show the current tab, and add an "active" class to the button that opened the tab
    document.getElementById(tabName).style.display = "block";
    evt.currentTarget.className += " tablinks-active";
} 


// used to match inputted text against list of gene names
function autocompleteMatch(input) {
    if (input == '') {
        return [];
    }
    var reg = new RegExp(input, "i") // "i" flag is for ignoring letter-case in gene search
    return geneList.filter(function(term) {
        if (term.match(reg)) {
        return term;
        }
    });
}

// displays genes matching regex for search term
function showResults(val) {
    //console.log(`showResults(${val})`);
    let res = document.getElementById("genesearchresult");
    res.innerHTML = '';
    let list = '';
    let terms = autocompleteMatch(val);
    terms.sort((a, b) => a.length - b.length); // sorting without a custom function will sort alphabetically
    //console.log('terms.length: ', terms.length, 'terms: ', terms);
    if(terms.length > 5){// set a max of 5 listed genes
        for (let i=0; i<5; i++) { 
            list += `<li class="dropdown-item" onclick="toogleInHighlightList('${terms[i]}')">` + terms[i] + '</li>';
        }
        list += `<li class="dropdown-item">` + '...' + '</li>';
    }
    else{
        for (const element of terms) { 
            list += `<li class="dropdown-item" onclick="toogleInHighlightList('${element}')">` + element + '</li>';
        }
    }
    res.innerHTML = '<ul class="list-group">' + list + '</ul>';
}

// updates the list of genes to highlight, as well as the dom to display these genes
function toogleInHighlightList(val, listOptions, listSelected){
    //console.log(`toogleInHighlightList(${val})`);
    // update gene list to highlight
    if(geneListToHighlight.includes(val)){
        geneListToHighlight = geneListToHighlight.filter( function(v) {
            return v !== val
        });
    }
    else{
        geneListToHighlight.push(val);
    }
    //console.log('new gene list to highlight: ', geneListToHighlight);

    // update UI
    let res = document.getElementById("genehighlightlist");
    res.innerHTML = '';
    let list = '';
    for (let i=0; i<geneListToHighlight.length; i++) {
        list += `<li class="dropdown-item">` + geneListToHighlight[i] + '</li>';
    }
    res.innerHTML = '<ul class="list-group" id="highlightlistford3">' + list + '</ul>';

    // update d3

}
