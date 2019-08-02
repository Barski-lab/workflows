cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_root = function(basename) {
          return basename.split('.').slice(0,1).join('.');
      };


inputs:

  fastq_file_1:
    type: File
    doc: "Paired-end sequencing data 1 in FASTQ format (fastq, fq, bzip2, gzip, zip)"

  fastq_file_2:
    type: File
    doc: "Paired-end sequencing data 2 in FASTQ format (fastq, fq, bzip2, gzip, zip)"

  indices_folder:
    type: Directory
    doc: "Directory with the genome indices generated by Bowtie2"

  exclude_chromosome:
    type: string?
    default: "chrM chrY chrX"
    doc: "Case-sensitive space-separated chromosome list to be excluded"

  blacklisted_regions_bed:
    type: File
    doc: "Blacklisted genomic regions file in BED format"

  chrom_length_file:
    type: File
    doc: "Chromosome length file"

  genome_size:
    type: string
    doc: "The length of the mappable genome (hs, mm, ce, dm or number, for example 2.7e9)"

  genome_fasta_file:
    type: File
    secondaryFiles: $(self.basename+".fai")  # due to bug in cwltool==1.0.20190621234233
    doc: "Reference genome sequence FASTA and FAI index files"

  threads:
    type: int?
    default: 4
    doc: "Number of threads for those steps that support multithreading"


outputs:

  fastq_1_qc_report_original:
    type: File
    outputSource: rename_qc_fastq_1_report/target_file

  fastq_1_adapter_trimming_report:
    type: File
    outputSource: trim_adapters/report_file

  fastq_2_qc_report_original:
    type: File
    outputSource: rename_qc_fastq_2_report/target_file

  fastq_2_adapter_trimming_report:
    type: File
    outputSource: trim_adapters/report_file_pair
  
  alignment_log:
    type: File
    outputSource: align_reads/output_log

  aligned_reads:
    type: File
    outputSource: sort_and_index/bam_bai_pair

  alignment_statistics:
    type: File
    outputSource: get_alignment_statistics/log_file

  read_redundancy_estimation:
    type: File
    outputSource: estimate_read_redundancy/estimates_file

  genome_coverage:
    type: File
    outputSource: convert_genome_coverage_to_bigwig/bigwig_file

  peak_calling_log:
    type: File
    outputSource: call_peaks/macs_log

  merged_peaks_with_counts:
    type: File
    outputSource: count_tags/intersected_file

  merged_peaks_sequences:
    type: File
    outputSource: get_sequences/sequences_file


steps:


# -----------------------------------------------------------------------------------


  extract_fastq_1:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_1
    out: [fastq_file]

  qc_fastq_1:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_fastq_1/fastq_file
    out:
      - summary_file
      - html_file

  rename_qc_fastq_1_report:
    run: ../tools/rename.cwl
    in:
      source_file: qc_fastq_1/html_file
      target_filename:
        source: fastq_file_1
        valueFrom: $(get_root(self.basename)+"_qc_original.html")
    out: [target_file]

  trigger_fastq_1_adapter_trimming:
    run: ../expressiontools/fastqc-results-trigger.cwl
    in:
      summary_file: qc_fastq_1/summary_file
    out: [trigger]


# -----------------------------------------------------------------------------------


  extract_fastq_2:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_2
    out: [fastq_file]

  qc_fastq_2:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_fastq_2/fastq_file
    out:
      - summary_file
      - html_file

  rename_qc_fastq_2_report:
    run: ../tools/rename.cwl
    in:
      source_file: qc_fastq_2/html_file
      target_filename:
        source: fastq_file_2
        valueFrom: $(get_root(self.basename)+"_qc_original.html")
    out: [target_file]

  trigger_fastq_2_adapter_trimming:
    run: ../expressiontools/fastqc-results-trigger.cwl
    in:
      summary_file: qc_fastq_2/summary_file
    out: [trigger]


# -----------------------------------------------------------------------------------


  trim_adapters:
    run: ../tools/trimgalore.cwl
    in:
      trigger:
        source: [trigger_fastq_1_adapter_trimming/trigger, trigger_fastq_2_adapter_trimming/trigger]
        valueFrom: $(self[0] || self[1])               # run trimgalore if at least one of input fastq files failed quality check
      input_file: extract_fastq_1/fastq_file
      input_file_pair: extract_fastq_2/fastq_file
      quality:
        default: 30      # Why do we need it if default should be 20
      dont_gzip:
        default: true    # should make it faster
      length:
        default: 30      # discard all reads shorter than 30 bp
      paired:
        default: true
    out:
      - trimmed_file
      - trimmed_file_pair
      - report_file
      - report_file_pair


# -----------------------------------------------------------------------------------


  rename_trimmed_fastq_1:
    run: ../tools/rename.cwl
    in:
      source_file: trim_adapters/trimmed_file
      target_filename:
        source: fastq_file_1
        valueFrom: $(get_root(self.basename) + ".fastq")
    out: [target_file]

  rename_trimmed_fastq_2:
    run: ../tools/rename.cwl
    in:
      source_file: trim_adapters/trimmed_file_pair
      target_filename:
        source: fastq_file_2
        valueFrom: $(get_root(self.basename) + ".fastq")
    out: [target_file]


# -----------------------------------------------------------------------------------


  align_reads:
    run: ../tools/bowtie2.cwl
    in:
      filelist: rename_trimmed_fastq_1/target_file
      filelist_mates: rename_trimmed_fastq_2/target_file
      indices_folder: indices_folder
      end_to_end_very_sensitive:
        default: true
      maxins:
        default: 2000
      no_discordant:          # do we need it?
        default: true
      no_mixed:               # do we need it?
        default: true
      threads: threads
    out:
      - output
      - output_log

  sort_and_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: align_reads/output
      threads: threads
    out: [bam_bai_pair]

  get_alignment_statistics:
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: sort_and_index/bam_bai_pair
    out:
      - log_file
      - average_length
      - reads_mapped

  estimate_read_redundancy:
    run: ../tools/preseq-lc-extrap.cwl
    in:
      bam_file: sort_and_index/bam_bai_pair
      pe_mode:
        default: true
      extrapolation:
        default: 1000000000
    out: [estimates_file]


# -----------------------------------------------------------------------------------


  filter_reads:
    run: ../tools/samtools-filter.cwl
    in:
      bam_bai_pair: sort_and_index/bam_bai_pair
      exclude_chromosome: exclude_chromosome
      quality:
        default: 30                 # how do we define 30 (range is from 0 to 255)
    out: [filtered_bam_bai_pair]
  
  remove_duplicates:
    run: ../tools/samtools-rmdup.cwl
    in:
      bam_file: filter_reads/filtered_bam_bai_pair
    out: [rmdup_output]

  convert_bam_to_bed:
    run: ../tools/bedtools-bamtobed.cwl
    in:
      bam_file: remove_duplicates/rmdup_output   # do we need to split reads by N
    out: [bed_file]

  shift_reads:
    run: ../tools/custom-bash.cwl
    in:
      input_file: convert_bam_to_bed/bed_file
      script:
        default: cat "$0" | awk 'BEGIN {OFS = "\t"}; {if ($6 == "+") print $1,$2+4,$3+4,$4,$5,$6; else print $1,$2-5,$3-5,$4,$5,$6}' > `basename $0`
    out: [output_file]

  remove_blacklisted:
    run: ../tools/bedtools-intersect.cwl
    in:
      file_a: shift_reads/output_file
      file_b: blacklisted_regions_bed
      no_overlaps:
        default: true
    out: [intersected_file]
                

# -----------------------------------------------------------------------------------


  group_by_chromosome:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: remove_blacklisted/intersected_file
      key:
        default: ["1,1"]
    out: [sorted_file]

  get_genome_coverage:
    run: ../tools/bedtools-genomecov.cwl
    in:
      input_file: group_by_chromosome/sorted_file
      chrom_length_file: chrom_length_file
      depth:
        default: "-bg"
      mapped_reads_number: get_alignment_statistics/reads_mapped
    out: [genome_coverage_file]

  sort_genome_coverage:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: get_genome_coverage/genome_coverage_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted_file]

  convert_genome_coverage_to_bigwig:
    run: ../tools/ucsc-bedgraphtobigwig.cwl
    in:
      bedgraph_file: sort_genome_coverage/sorted_file
      chrom_length_file: chrom_length_file
    out: [bigwig_file]


# -----------------------------------------------------------------------------------


  call_peaks:
    run: ../tools/macs2-callpeak.cwl
    in:
      treatment_file: remove_blacklisted/intersected_file
      format_mode:
        default: "BED"
      genome_size: genome_size
      keep_dup:
        default: "all"
      nomodel:
        default: true
      shift:
        source: get_alignment_statistics/average_length
        valueFrom: $(-Math.round(self/2))
      extsize: get_alignment_statistics/average_length
    out:
      - narrow_peak_file
      - macs_log

  sort_peaks:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: call_peaks/narrow_peak_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted_file]

  merge_peaks:
    run: ../tools/bedtools-merge.cwl
    in:
      bed_file: sort_peaks/sorted_file
    out: [merged_bed_file]

  count_tags:
    run: ../tools/bedtools-intersect.cwl
    in:
      file_a: merge_peaks/merged_bed_file
      file_b: remove_blacklisted/intersected_file
      count:
        default: true
    out: [intersected_file]

  get_sequences:
    run: ../tools/bedtools-getfasta.cwl
    in:
      genome_fasta_file: genome_fasta_file
      intervals_file: merge_peaks/merged_bed_file
    out: [sequences_file]