cwlVersion: v1.0
class: Workflow


inputs:

  fastq_file:
    type: File
    format: "http://edamontology.org/format_1930"
    label: "FASTQ file"
    doc: "Uncompressed or gzipped FASTQ file, single-end"

  adapters_file:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Adapters file"
    doc: "FASTA file containing adapters"

  indices_folder:
    type: Directory
    label: "Bismark indices folder"
    doc: "Path to Bismark generated indices folder"

  processes:
    type: int?
    label: "Number of Bismark instances to run"
    doc: |
      "Set the number of parallel Bismark instances to run concurrently.
       Each Bismark instance simultainously runs the methylation extractor,
       samtools stream and GZIP streams"

  threads:
    type: int?
    label: "Number of Bowtie2/Trimmomatic/Samtools threads to use"
    doc: "Set the number of threads for Bowtie2, Trimmomatic, Samtools"


outputs:

  bambai_pair:
    type: File
    label: "BAM alignment and BAI index files"
    doc: "Bismark generated coordinate sorted BAM alignment and BAI index files"
    outputSource: samtools_sort_index/bam_bai_pair

  bismark_alignment_report:
    type: File
    label: "Bismark alignment and methylation report"
    doc: "Bismark generated alignment and methylation summary report"
    outputSource: bismark_align/alignment_report

  chg_context_file:
    type: File
    label: "CHG methylation call"
    doc: "CHG methylation call"
    outputSource: bismark_extract_methylation/chg_context_file

  chh_context_file:
    type: File
    label: "CHH methylation call"
    doc: "CHH methylation call"
    outputSource: bismark_extract_methylation/chh_context_file

  cpg_context_file:
    type: File
    label: "CpG methylation call"
    doc: "CpG methylation call"
    outputSource: bismark_extract_methylation/cpg_context_file

  mbias_plot:
    type: File
    label: "Methylation bias plot"
    doc: "QC data showing methylation bias across read lengths"
    outputSource: bismark_extract_methylation/mbias_plot

  mbias_plot_png:
    type: File
    label: "Methylation bias plot (PNG)"
    doc: "QC data showing methylation bias across read lengths"
    outputSource: bismark_extract_methylation/mbias_plot_png

  bedgraph_coverage_file:
    type: File
    label: "Methylation statuses bedGraph coverage file"
    doc: "Coverage text file summarising cytosine methylation values in bedGraph format (tab-delimited; 0-based start coords, 1-based end coords)"
    outputSource: bismark_extract_methylation/bedgraph_coverage_file

  bismark_coverage_file:
    type: File
    label: "Methylation statuses Bismark coverage file"
    doc: "Coverage text file summarising cytosine methylation values in Bismark format (tab-delimited, 1-based genomic coords)"
    outputSource: bismark_extract_methylation/bismark_coverage_file

  genome_wide_methylation_report:
    type: File
    label: "Genome-wide cytosine methylation report"
    doc: "Genome-wide methylation report for all cytosines in the genome"
    outputSource: bismark_extract_methylation/genome_wide_methylation_report

  splitting_report:
    type: File
    label: "Methylation extraction log"
    doc: "Log file giving summary statistics about methylation extraction"
    outputSource: bismark_extract_methylation/splitting_report

  collected_report:
    type: File
    label: "HTML report page"
    doc: "Bismark generated graphical HTML report page"
    outputSource: bismark_report/collected_report


steps:

  trim_adapters:
    run: ../tools/trimmomatic.cwl
    in:
      fastq_file_upstream: fastq_file
      adapters_file: adapters_file
      lib_type:
        default: "SE"
      illuminaclip_step_param:
        default: "2:30:10"
      leading_step:
        default: 30
      trailing_step:
        default: 30
      sliding_window_step:
        default: "5:30"
      minlen_step:
        default: 25
      threads: threads
    out: [upstream_trimmed_file]

  bismark_align:
    run: ../tools/bismark-align.cwl
    in:
      fastq_file: trim_adapters/upstream_trimmed_file
      indices_folder: indices_folder
      processes: processes
      threads: threads
    out: [bam_file, alignment_report]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: bismark_align/bam_file
      threads: threads
    out: [bam_bai_pair]

  bismark_extract_methylation:
    run: ../tools/bismark-extract-methylation.cwl
    in:
      genome_folder: indices_folder
      bam_file: samtools_sort_index/bam_bai_pair
      processes: processes
      threads: threads
    out:
      - chg_context_file
      - chh_context_file
      - cpg_context_file
      - mbias_plot
      - mbias_plot_png
      - bedgraph_coverage_file
      - bismark_coverage_file
      - genome_wide_methylation_report
      - splitting_report

  bismark_report:
    run: ../tools/bismark-report.cwl
    in:
      alignment_report: bismark_align/alignment_report
      splitting_report: bismark_extract_methylation/splitting_report
      mbias_report: bismark_extract_methylation/mbias_plot
    out: [collected_report]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "bismark-methylation-se"
s:downloadUrl: https://github.com/Barski-lab/workflows/blob/master/workflows/bismark-methylation-se.cwl
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
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Bismark Methylation pipeline. We can use indices_folder as genome_folder for bismark_extract_methylation step,
  because it insludes the original FASTA files too.

s:about: |

  Sequence reads are first transformed into fully bisulfite-converted forward (C->T) and reverse read (G->A conversion of
  the forward strand) versions, before they are aligned to similarly converted versions of the genome (also C->T and G->A
  converted). Sequence reads that produce a unique best alignment from the four alignment processes against the bisulfite
  genomes (which are running in parallel) are then compared to the normal genomic sequence and the methylation state of
  all cytosine positions in the read is inferred. A read is considered to align uniquely if an alignment has a unique best
  alignment score (as reported by the AS:i field). If a read produces several alignments with the same number of mismatches
  or with the same alignment score (AS:i field), a read (or a read-pair) is discarded altogether.





