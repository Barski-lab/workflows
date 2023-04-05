cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_label = function(i) {
        var rootname = inputs.barcode_metrics_report[i].basename.split('.').slice(0,-1).join('.');
        rootname = (rootname=="")?inputs.barcode_metrics_report[i].basename:rootname;
        return inputs.gem_well_labels?inputs.gem_well_labels[i].replace(/\t|\s|\[|\]|\>|\<|,|\./g, "_"):rootname;
    };
- class: InitialWorkDirRequirement
  listing: |
    ${
      var entry = "library_id,fragments,cells\n"
      for (var i=0; i < inputs.barcode_metrics_report.length; i++){
        entry += get_label(i) + "," + inputs.fragments_file_from_count[i].path + "," + inputs.barcode_metrics_report[i].path + "\n";
      }
      return [{
        "entry": entry,
        "entryname": "metadata.csv"
      }];
    }


hints:
- class: DockerRequirement
  dockerPull: cumulusprod/cellranger-atac:2.1.0


inputs:

  fragments_file_from_count:
    type: File[]
    secondaryFiles:
    - .tbi
    doc: |
      Array of files containing count and barcode information for
      every ATAC fragment observed in the "cellranger-atac count"
      experiment in TSV format.

  barcode_metrics_report:
    type: File[]
    doc: |
      Array of files with per-barcode fragment counts & metrics
      produced by "cellranger-atac count" command in CSV format

  gem_well_labels:
    type:
    - "null"
    - string[]
    doc: |
      Array of GEM well identifiers to be used for labeling purposes only.
      If not provided use rootnames of files from the barcode_metrics_report
      input

  indices_folder:
    type: Directory
    inputBinding:
      position: 5
      prefix: "--reference"
    doc: |
      Path to folder containing a Cell Ranger ATAC or Cell Ranger
      ARC reference. Should be generated by "cellranger-atac mkref"
      or "cellranger-arc mkref" commands

  normalization_mode:
    type:
    - "null"
    - type: enum
      name: "normalization"
      symbols: ["none", "depth"]
    inputBinding:
      position: 6
      prefix: "--normalize"
    doc: |
      Library depth normalization mode: depth, none.
      Default: depth

  threads:
    type: int?
    inputBinding:
      position: 7
      prefix: "--localcores"
    doc: |
      Set max cores the pipeline may request at one time.
      Default: all available

  memory_limit:
    type: int?
    inputBinding:
      position: 8
      prefix: "--localmem"
    doc: |
      Set max GB the pipeline may request at one time
      Default: all available

  virt_memory_limit:
    type: int?
    inputBinding:
      position: 9
      prefix: "--localvmem"
    doc: |
      Set max virtual address space in GB for the pipeline
      Default: all available


outputs:

  web_summary_report:
    type: File
    outputBinding:
      glob: "aggregated/outs/web_summary.html"
    doc: |
      Run summary metrics and charts in HTML format

  metrics_summary_report_json:
    type: File
    outputBinding:
      glob: "aggregated/outs/summary.json"
    doc: |
      Run summary metrics in JSON format

  metrics_summary_report_csv:
    type: File
    outputBinding:
      glob: "aggregated/outs/summary.csv"
    doc: |
      Run summary metrics in CSV format

  barcode_metrics_report:
    type: File
    outputBinding:
      glob: "aggregated/outs/singlecell.csv"
    doc: |
      Per-barcode fragment counts & metrics in CSV format

  fragments_file:
    type: File
    outputBinding:
      glob: "aggregated/outs/fragments.tsv.gz"
    secondaryFiles:
    - .tbi
    doc: |
      Count and barcode information for every ATAC fragment observed
      in the aggregated experiment in TSV format

  peaks_bed_file:
    type: File
    outputBinding:
      glob: "aggregated/outs/peaks.bed"
    doc: |
      Locations of open-chromatin regions identified in the
      aggregated experiment (these regions are referred to
      as "peaks")

  peak_annotation_file:
    type: File
    outputBinding:
      glob: "aggregated/outs/peak_annotation.tsv"
    doc: |
      Annotations of peaks based on genomic proximity alone

  secondary_analysis_report_folder:
    type: Directory
    outputBinding:
      glob: "aggregated/outs/analysis"
    doc: |
      Folder with secondary analysis results

  filtered_feature_bc_matrix_folder:
    type: Directory
    outputBinding:
      glob: "aggregated/outs/filtered_peak_bc_matrix"
    doc: |
      Folder with aggregated filtered peak-barcode matrices
      containing only cellular barcodes in MEX format.

  filtered_feature_bc_matrix_h5:
    type: File
    outputBinding:
      glob: "aggregated/outs/filtered_peak_bc_matrix.h5"
    doc: |
      Aggregated filtered peak-barcode matrices containing
      only cellular barcodes in HDF5 format.

  filtered_tf_bc_matrix_folder:
    type: Directory
    outputBinding:
      glob: "aggregated/outs/filtered_tf_bc_matrix"
    doc: |
      Folder with aggregated filtered tf-barcode matrices
      containing only cellular barcodes in MEX format.

  filtered_tf_bc_matrix_h5:
    type: File
    outputBinding:
      glob: "aggregated/outs/filtered_tf_bc_matrix.h5"
    doc: |
      Aggregated filtered tf-barcode matrices containing
      only cellular barcodes in HDF5 format.

  aggregation_metadata:
    type: File
    outputBinding:
      glob: "aggregated/outs/aggregation_csv.csv"
    doc: |
      Aggregation CSV file

  loupe_browser_track:
    type: File
    outputBinding:
      glob: "aggregated/outs/cloupe.cloupe"
    doc: |
      Loupe Browser visualization and analysis file

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["cellranger-atac", "aggr", "--disable-ui", "--id", "aggregated", "--csv", "metadata.csv"]


stdout: cellranger_atac_aggr_stdout.log
stderr: cellranger_atac_aggr_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cellranger ATAC Aggregate"
s:name: "Cellranger ATAC Aggregate"
s:alternateName: "Aggregates outputs from multiple runs of Cell Ranger Count Chromatin Accessibility experiments"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/cellranger-atac-aggr.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Cellranger ATAC Aggregate

  Aggregates outputs from multiple runs of Cell Ranger Count Chromatin
  Accessibility experiments

  Parameters set by default:
  --disable-ui - no need in any UI when running in Docker container
  --id         - hardcoded to `aggregated` as we want to return the
                 content of the outputs folder as separate outputs

  Skipped parameters:
  --description
  --peaks
  --nosecondary
  --dim-reduce
  --dry
  --jobmode
  --mempercore
  --maxjobs
  --jobinterval
  --overrides
  --uiport
  --noexit
  --nopreflight


s:about: |
  USAGE:
    cellranger-atac aggr [OPTIONS] --id <ID> --csv <CSV> --reference <PATH>

  OPTIONS:
    --id <ID>
        A unique run id and output folder name [a-zA-Z0-9_-]+ of maximum length 64 characters
    --description <TEXT>
        Sample description to embed in output files
        [default: ]
    --csv <CSV>
        Path to CSV file enumerating `cellranger-atac count` outputs.
        For example, a CSV for aggregating two samples would look as follows (blank lines are ignored):
        library_id,fragments,cells
        L1,/data/L1/outs/fragments.tsv.gz,/data/L1/outs/singlecell.csv
        L2,/data/L2/outs/fragments.tsv.gz,/data/L2/outs/singlecell.csv
        Optionally, metadata associated with these libraries can be specified using additional columns. This information is not used by the pipeline but will be available in
        the Loupe file for visualization.
    --reference <PATH>
        Path to folder containing a Cell Ranger ATAC or Cell Ranger ARC reference
    --peaks <BED>
        Override peak caller: specify peaks to use in downstream analyses from supplied 3-column BED file. The supplied peaks file must be sorted by position and not contain
        overlapping peaks; comment lines beginning with `#` are allowed
    --normalize <MODE>
        Library depth normalization mode
        [default: depth]
        [possible values: depth, none]
    --nosecondary
        Disable secondary analysis, e.g. clustering
    --dim-reduce <STR>
        Dimensionality reduction mode for clustering
        [default: lsa]
        [possible values: lsa, pca, plsa]
    --dry
        Do not execute the pipeline. Generate a pipeline invocation (.mro) file and stop
    --jobmode <MODE>
        Job manager to use. Valid options: local (default), sge, lsf, slurm or path to a .template file. Search for help on "Cluster Mode" at support.10xgenomics.com for more
        details on configuring the pipeline to use a compute cluster
        [default: local]
    --localcores <NUM>
        Set max cores the pipeline may request at one time. Only applies to local jobs
    --localmem <NUM>
        Set max GB the pipeline may request at one time. Only applies to local jobs
    --localvmem <NUM>
        Set max virtual address space in GB for the pipeline. Only applies to local jobs
    --mempercore <NUM>
        Reserve enough threads for each job to ensure enough memory will be available, assuming each core on your cluster has at least this much memory available. Only applies
        to cluster jobmodes
    --maxjobs <NUM>
        Set max jobs submitted to cluster at one time. Only applies to cluster jobmodes
    --jobinterval <NUM>
        Set delay between submitting jobs to cluster, in ms. Only applies to cluster jobmodes
    --overrides <PATH>
        The path to a JSON file that specifies stage-level overrides for cores and memory. Finer-grained than --localcores, --mempercore and --localmem. Consult
        https://support.10xgenomics.com/ for an example override file
    --uiport <PORT>
        Serve web UI at http://localhost:PORT
    --disable-ui
        Do not serve the web UI
    --noexit
        Keep web UI running after pipestance completes or fails
    --nopreflight
        Skip preflight checks
    --help
        Print help information