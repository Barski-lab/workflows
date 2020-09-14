cwlVersion: v1.0
class: Workflow


requirements:
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
    doc: "Directory with the genome indices generated by STAR"

  blacklisted_regions_bed:
    type: File
    doc: "Blacklisted genomic regions file in BED format"

  genome_fasta_file:
    type: File
    secondaryFiles:
    - .fai
    doc: "Reference genome sequence FASTA and FAI index files"

  genome_size:
    type: string
    doc: "The length of the mappable genome (hs, mm, ce, dm or number, for example 2.7e9)"    

  exclude_chromosome:
    type: string?
    default: "chrM chrY chrX"
    doc: "Case-sensitive space-separated chromosome list to be excluded"    

  trim_adapter_criteria:  # TODO: see what criteria we need to trigger adapter trimming
    type: string?
    default: ".*Per base sequence quality.*|.*Per sequence quality scores.*|.*Overrepresented sequences.*|.*Adapter Content.*"
    doc: "Regex for setting failed FastQC summary items any of which will trigger adapter trimming"

  threads:
    type: int?
    default: 4
    doc: "Number of threads for those steps that support multithreading"


outputs:

  fastqc_report_fastq_1:
    type: File
    outputSource: rename_fastqc_report_fastq_1/target_file

  fastqc_report_fastq_2:
    type: File
    outputSource: rename_fastqc_report_fastq_2/target_file

  adapter_trimming_report_fastq_1:
    type: File
    outputSource: rename_adapter_trimming_report_fastq_1/target_file

  adapter_trimming_report_fastq_2:
    type: File
    outputSource: rename_adapter_trimming_report_fastq_2/target_file
  
  star_alignment_statistics_report:
    type: File
    outputSource: align_reads/log_final

  star_alignment_progress_report:
    type: File
    outputSource: align_reads/log_progress

  uniquely_mapped_sorted_reads_as_bam:
    type: File
    outputSource: sort_and_index/bam_bai_pair

  uniquely_mapped_sorted_reads_as_bam_statistics_report:
    type: File
    outputSource: get_bam_statistics/log_file

  samtools_deduplication_report:
    type: File
    outputSource: remove_duplicates/markdup_report

  uniquely_mapped_sorted_filtered_deduped_reads_as_bam:
    type: File
    outputSource: remove_duplicates/deduplicated_bam_bai_pair

  uniquely_mapped_sorted_filtered_deduped_reads_as_bam_statistics_report:
    type: File
    outputSource: get_bam_statistics_after_filtering/log_file

  filtered_deduplicated_reads_as_bed:
    type: File
    outputSource: convert_bam_to_bed/bed_file

  filtered_deduplicated_shifted_reads_as_bed:
    type: File
    outputSource: shift_reads/output_file

  filtered_deduplicated_sorted_shifted_reads_wo_blacklisted_as_bed:
    type: File
    outputSource: sort_bed/sorted_file

  genome_coverage_as_bigwig:
    type: File
    outputSource: convert_bedgraph_to_bigwig/bigwig_file

  macs2_peak_calling_report:
    type: File
    outputSource: call_peaks/macs_log

  macs2_called_peaks_as_sorted_narrow_peak:
    type: File
    outputSource: sort_peaks/sorted_file

  macs2_called_peaks_merged_sorted_narrow_peak:
    type: File
    outputSource: sort_merged_peaks/sorted_file

  tag_counts_within_merged_sorted_peaks:
    type: File
    outputSource: count_tags/intersected_file

  sequences_within_merged_sorted_peaks:
    type: File
    outputSource: get_sequences/sequences_file


steps:

# -----------------------------------------------------------------------------------

  extract_fastq_1:
    doc: |
      Tries to uncompress input fastq file 1. If file wasn't
      compressed, returns it unchanged
    run: ../../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_1
    out:
    - fastq_file                                                      # Uncompressed fastq file 1

  run_fastqc_fastq_1:
    doc: Runs FastQC for uncompressed fastq file 1
    run: ../../tools/fastqc.cwl
    in:
      reads_file: extract_fastq_1/fastq_file
    out:
    - summary_file                                                    # Will be used to check if adapter trimming needed
    - html_file                                                       # Will be renamed and returned as workflow output

  rename_fastqc_report_fastq_1:
    doc: |
      Renames FastQC report for fastq file 1 (technical step)
    run: ../../tools/rename.cwl
    in:
      source_file: run_fastqc_fastq_1/html_file
      target_filename:
        source: fastq_file_1
        valueFrom: $(get_root(self.basename)+"_fastqc_report.html")   # Use root name of the original fastq file 1 and updated suffix
    out:
    - target_file                                                     # Returned as workflow output. Not used in any other calculations

  trigger_adapter_trimming_fastq_1:
    doc: |
      Checks if we want to trim adaprters based on FastQC report
      for fastq file 1
    run: ../../tools/fastqc-results-trigger.cwl
    in:
      summary_file: run_fastqc_fastq_1/summary_file
      criteria: trim_adapter_criteria
    out:
    - trigger                                                         # True, if at least one of the items from trim_adapter_criteria failed

# -----------------------------------------------------------------------------------

  extract_fastq_2:
    doc: |
      Tries to uncompress input fastq file 2. If file wasn't
      compressed, returns it unchanged
    run: ../../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_2
    out:
    - fastq_file                                                      # Uncompressed fastq file 2

  run_fastqc_fastq_2:
    doc: Runs FastQC for uncompressed fastq file 2
    run: ../../tools/fastqc.cwl
    in:
      reads_file: extract_fastq_2/fastq_file
    out:
    - summary_file                                                    # Will be used to check if adapter trimming needed
    - html_file                                                       # Will be renamed and returned as workflow output

  rename_fastqc_report_fastq_2:
    doc: Checks if we want to trim adaprters based on FastQC report for fastq file 2
    run: ../../tools/rename.cwl
    in:
      source_file: run_fastqc_fastq_2/html_file
      target_filename:
        source: fastq_file_2
        valueFrom: $(get_root(self.basename)+"_fastqc_report.html")   # Use root name of the original fastq file 2 and updated suffix
    out:
    - target_file                                                     # Returned as workflow output. Not used in any other calculations
  
  trigger_adapter_trimming_fastq_2:
    doc: Checks if we want to trim adaprters based on FastQC report for fastq file 2
    run: ../../tools/fastqc-results-trigger.cwl
    in:
      summary_file: run_fastqc_fastq_2/summary_file
      criteria: trim_adapter_criteria
    out:
    - trigger                                                         # True, if at least one of the items from trim_adapter_criteria failed

# -----------------------------------------------------------------------------------

  trim_adapters:
    doc: |
      Trims adapters when at least on of the input fastq files
      failed quality control
    run: ../../tools/trimgalore.cwl
    in:
      trigger:
        source:
        - trigger_adapter_trimming_fastq_1/trigger
        - trigger_adapter_trimming_fastq_2/trigger
        valueFrom: $(self[0] || self[1])                              # OR 
      input_file: extract_fastq_1/fastq_file
      input_file_pair: extract_fastq_2/fastq_file
      quality:
        default: 30                                                   # TODO: Do we really need 30? The default was 20
      dont_gzip:
        default: true                                                 # Don't need to compress (makes workflow run faster)
      length:
        default: 30                                                   # Discards both reads if any of them became shorter than 30 bp after adapter trimming
      paired:
        default: true                                                 # Tells TrimGalore that we work with paired-end data
    out:
    - trimmed_file                                                    # Trimmed uncompressed fastq file 1
    - trimmed_file_pair                                               # Trimmed uncompressed fastq file 2
    - report_file                                                     # Renamed and returned as workflow output. Not used in any other calculations
    - report_file_pair                                                # Renamed and returned as workflow output. Not used in any other calculations

  rename_adapter_trimming_report_fastq_1:
    doc: |
      Renames adapter trimming report for fastq file 1
      (technical step)
    run: ../../tools/rename.cwl
    in:
      source_file: trim_adapters/report_file                          # Adapter trimming report for fastq file 1
      target_filename:                                                # Use basename of the original fastq file 1
        source: fastq_file_1
        valueFrom: $(get_root(self.basename) + "_adapter_trimming_report.txt")
    out:
    - target_file                                                     # Renamed adapter trimming report for fastq file 1

  rename_adapter_trimming_report_fastq_2:
    doc: |
      Renames adapter trimming report for fastq file 2
      (technical step)
    run: ../../tools/rename.cwl
    in:
      source_file: trim_adapters/report_file_pair                     # Adapter trimming report for fastq file 2
      target_filename:                                                # Use basename of the original fastq file 2
        source: fastq_file_2
        valueFrom: $(get_root(self.basename) + "_adapter_trimming_report.txt")
    out:
    - target_file                                                     # Renamed adapter trimming report for fastq file 2


  rename_trimmed_fastq_1:
    doc: |
      Renames trimmed fastq file 1 (technical step to get rid of
      the suffix added by TrimGalore to the filename)
    run: ../../tools/rename.cwl
    in:
      source_file: trim_adapters/trimmed_file                         # Trimmed uncompressed fastq file 1
      target_filename:
        source: fastq_file_1
        valueFrom: $(get_root(self.basename) + ".fastq")              # Use basename of the original fastq file 1
    out:
    - target_file                                                     # Renamed trimmed uncompressed fastq file 1

  rename_trimmed_fastq_2:
    doc: |
      Renames trimmed fastq file 2 (technical step to get rid of
      the suffix added by TrimGalore to the filename)
    run: ../../tools/rename.cwl
    in:
      source_file: trim_adapters/trimmed_file_pair                    # Trimmed uncompressed fastq file 2
      target_filename:
        source: fastq_file_2
        valueFrom: $(get_root(self.basename) + ".fastq")              # Use basename of the original fastq file 2
    out:
    - target_file                                                     # Renamed trimmed uncompressed fastq file 2

# -----------------------------------------------------------------------------------

  align_reads:
    doc: |
      Aligns trimmed uncompressed fastq files using STAR. Skips all
      multimapped reads. Unmapped reads are not reported.
    run: ../../tools/star-alignreads.cwl
    in:
      readFilesIn:
      - rename_trimmed_fastq_1/target_file
      - rename_trimmed_fastq_2/target_file
      genomeDir: indices_folder
      outSAMtype:
        default:
        - "SAM"                                                       # Returns unsorted reads in SAM format
      outSAMunmapped:
        default: "None"                                               # Do not report unmapped reads
      alignIntronMax:                                                 # TODO: Is it something specific to make STAR work for ATAC-Seq?
        default: 1
      alignEndsType:                                                  # TODO: Is it something specific to make STAR work for ATAC-Seq?
        default: "EndToEnd"
      alignMatesGapMax:                                               # TODO: Is it something specific to make STAR work for ATAC-Seq?
        default: 2000        
      outFilterMultimapNmax:
        default: 1                                                    # Do not report multimapped reads
      outFilterMismatchNmax:
        default: 5                                                    # Allow maximum 5 mismatches per read pair
      outFileNamePrefix:                                              # Sets the output files prefix to be a combination of fastq files 1,2 root names
        source:
        - rename_trimmed_fastq_1/target_file
        - rename_trimmed_fastq_2/target_file
        valueFrom: $(get_root(self[0].basename) + "_" + get_root(self[1].basename) + ".")
      threads: threads
    out:
    - aligned_file                                                    # Unsorted SAM file with ONLY uniquely mapped reads with max 5 mismatches per read pair
    - log_final                                                       # Returned as workflow output. Not used in any other calculations
    - log_progress                                                    # Returned as workflow output. Not used in any other calculations

  sort_and_index:
    doc: |
      Sorts and indexes SAM file. Returns coordinate sorted
      BAM file and BAI index
    run: ../../tools/samtools-sort-index.cwl
    in:
      sort_input: align_reads/aligned_file                            # Unsorted SAM file with ONLY uniquely mapped reads with max 5 mismatches per read pair
      sort_output_filename:                                           # Sets the output file name to be a combination of fastq files 1,2 root names
        source:
        - rename_trimmed_fastq_1/target_file
        - rename_trimmed_fastq_2/target_file
        valueFrom: $(get_root(self[0].basename) + "_" + get_root(self[1].basename) + "_uniquely_mapped_sorted.bam")
      threads: threads
    out:
    - bam_bai_pair                                                    # Coordinate sorted uniquely mapped reads with max 5 mismatches per read pair as BAM file and BAI index

  get_bam_statistics:
    doc: |
      Collects BAM statistics before any filters applied
    run: ../../tools/samtools-stats.cwl
    in:
      bambai_pair: sort_and_index/bam_bai_pair
      output_filename:                                                # Sets the output file name to include BAM file root name and updated suffix
        source: sort_and_index/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report.txt")
    out:
    - log_file                                                        # Returned as workflow output. Not used in any other calculations

# -----------------------------------------------------------------------------------

  filter_reads:
    doc: |
      Removes chromosomes we don't need. Returns coordinate sorted
      BAM file and BAI index
    run: ../../tools/samtools-filter.cwl
    in:
      bam_bai_pair: sort_and_index/bam_bai_pair                       # Coordinate sorted uniquely mapped reads with max 5 mismatches per read pair as BAM file and BAI index
      exclude_chromosome: exclude_chromosome                          # Space-separated chromosome list to be removed
      output_filename:                                                # Sets the output file name to include BAM file root name and updated suffix
        source: sort_and_index/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_filtered.bam")
    out:
    - filtered_bam_bai_pair                                           # Coordinate sorted filtered (with chromosomes removed) BAM file and BAI index

  get_bam_statistics_after_chrom_removal:
    doc: |
      Collects BAM statistics after unused chromosomes removal.
      We need it to get reads_mapped for genome coverage
      normalization
    run: ../../tools/samtools-stats.cwl
    in:
      bambai_pair: filter_reads/filtered_bam_bai_pair                 # Coordinate sorted filtered (with chromosomes removed) BAM file and BAI index
    out:
    - reads_mapped                                                    # mapped reads number will be used for genome coverage normaization

  remove_duplicates:
    doc: |
      Removes PCR duplicates. Returns coordinate sorted BAM file
      and BAI index
    run: ../../tools/samtools-markdup.cwl
    in:
      bam_bai_pair: filter_reads/filtered_bam_bai_pair                # Coordinate sorted filtered (with chromosomes removed) BAM file and BAI index
      threads: threads
      output_filename:                                                # Sets the output file name to include BAM file root name and updated suffix
        source: filter_reads/filtered_bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_deduped.bam")      
    out:
    - deduplicated_bam_bai_pair                                       # Coordinate sorted filtered (with chromosomes removed) and deduplicated BAM file and BAI index
    - markdup_report                                                  # Returned as workflow output. Not used in any other calculations

  get_bam_statistics_after_filtering:
    doc: |
      Collects BAM statistics after all filtering steps applied
    run: ../../tools/samtools-stats.cwl
    in:
      bambai_pair: remove_duplicates/deduplicated_bam_bai_pair        # Coordinate sorted filtered (with chromosomes removed) and deduplicated BAM file and BAI index
      output_filename:
        source: remove_duplicates/deduplicated_bam_bai_pair           # Sets the output file name to include BAM file root name and updated suffix
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report.txt")
    out:
    - log_file                                                        # Returned as workflow output. Not used in any other calculations
    - average_length                                                  # Used by MACS2 in peak calling

  convert_bam_to_bed:
    doc: |
      Converts coordinate sorted filtered deduplicated BAM to
      BED format. All paired-ends become single-reads (regions).
      Reads sequences are removed. We don't split reads by N or D
      in CIGAR (do we really have them in ATAC?)
    run: ../../tools/bedtools-bamtobed.cwl
    in:
      bam_file: remove_duplicates/deduplicated_bam_bai_pair           # Coordinate sorted filtered (with chromosomes removed) and deduplicated BAM file and BAI index
    out:
    - bed_file                                                        # Uniquely mapped filtered deduped reads as regions. Connections between pairs are lost

  shift_reads:
    doc: |
      Shifts and center uniquely mapped filtered deduped reads
      (regions) to Tn5 binding sites. Returns sorted by -k1,1
      -k2,2n -k3,3n results. Each region becomes 40bp long.
    run: ../../tools/custom-bash.cwl
    in:
      input_file: convert_bam_to_bed/bed_file                         # Uniquely mapped filtered deduped reads as regions
      param:
        source: convert_bam_to_bed/bed_file                           # Sets the output file name to include BED file root name and updated suffix
        valueFrom: $(get_root(self.basename)+"_shifted.bed")
      script:
        default: cat "$0" | awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, ($3 + 4) - 20, ($3 + 4) + 20, $4, $5, $6; else print $1, ($2 - 5) - 20, ($2 - 5) + 20, $4, $5, $6}' | sort -k1,1 -k2,2n -k3,3n> $1
    out:
    - output_file                                                     # Sorted by coordinates Tn5 binding sites

  remove_blacklisted:
    doc: |
      Removes all Tn5 binding sites that intersect blacklisted
      regions
    run: ../../tools/bedtools-intersect.cwl
    in:
      file_a: shift_reads/output_file                                 # Sorted by coordinates Tn5 binding sites
      file_b: blacklisted_regions_bed
      no_overlaps:                                                    # Reports those entries in file_a that have no overlaps with file_b
        default: true
    out:
    - intersected_file                                                # Tn5 binding sites without blacklisted regions

  sort_bed:
    doc: |
      Sorts Tn5 binding sites without blacklisted regions by -k1,1
      -k2,2n -k3,3n to be able to use it with bedtools genomecov
    run: ../../tools/linux-sort.cwl
    in:
      unsorted_file: remove_blacklisted/intersected_file              # Tn5 binding sites without blacklisted regions
      key:
        default: ["1,1","2,2n","3,3n"]
      output_filename:                                                # Sets the output file name to include BED file root name and updated suffix
        source: remove_blacklisted/intersected_file
        valueFrom: $(get_root(self.basename)+"_wo_blacklisted.bed")
    out:
    - sorted_file                                                     # Coordinate sorted Tn5 binding sites without blacklisted regions

  get_chr_name_length:
    doc: |
      Returns chrNameLength.txt file from the folder (technical step)
    run: ../../tools/get-file-by-name.cwl
    in:
      input_files: indices_folder                                     # STAR indices folder
      basename_regex:
        default: "chrNameLength.txt"
    out:
    - selected_file                                                   # chrNameLength.txt file

  convert_bed_to_bedgraph:
    doc: |
      Calculates genome coverage from Tn5 binding sites
      without blacklisted regions
    run: ../../tools/bedtools-genomecov.cwl
    in:
      input_file: sort_bed/sorted_file                                # Coordinate sorted Tn5 binding sites without blacklisted regions
      depth:
        default: "-bg"
      mapped_reads_number: get_bam_statistics_after_chrom_removal/reads_mapped  # Scaling coeficient. TODO: do we need to take mapped reads number before chromosomes removal?
      chrom_length_file: get_chr_name_length/selected_file
    out:
    - genome_coverage_file                                            # Genome coverage of Tn5 binding sites as bedGraph file

  sort_bedgraph:
    doc: |
      Returns coordinate sorted genome coverage of Tn5 binding sites
    run: ../../tools/linux-sort.cwl
    in:
      unsorted_file: convert_bed_to_bedgraph/genome_coverage_file     # Genome coverage of Tn5 binding sites as bedGraph file
      key:
        default: ["1,1","2,2n","3,3n"]
    out:
    - sorted_file                                                     # Coordinate sorted genome coverage of Tn5 binding sites as bedGraph file

  convert_bedgraph_to_bigwig:
    doc: |
      Converts coordinate sorted genome coverage of Tn5 binding
      sites to bigwig
    run: ../../tools/ucsc-bedgraphtobigwig.cwl
    in:
      bedgraph_file: sort_bedgraph/sorted_file                        # Coordinate sorted genome coverage of Tn5 binding sites as bedGraph file
      chrom_length_file: get_chr_name_length/selected_file            # chrNameLength.txt file
    out:
    - bigwig_file                                                     # Genome coverage of Tn5 binding sites as bigwig file

# -----------------------------------------------------------------------------------

  call_peaks:
    doc: |
      Calls peaks with MACS2. Need some explanation of
      used parameters
    run: ../../tools/macs2-callpeak.cwl
    in:
      treatment_file: sort_bed/sorted_file                            # Coordinate sorted genome coverage of Tn5 binding sites as bedGraph file
      format_mode:
        default: "BED"
      genome_size: genome_size
      keep_dup:
        default: "all"
      nomodel:
        default: true
      shift:                                                          # Use average read length from BAM file after all filters applied
        source: get_bam_statistics_after_filtering/average_length
        valueFrom: $(-Math.round(self/2))
      extsize: get_bam_statistics_after_filtering/average_length
    out:
    - narrow_peak_file                                                # Called peaks as not sorted narrowpeak file
    - macs_log                                                        # Returned as workflow output. Not used in any other calculations

  sort_peaks:
    doc: |
      Sorts called peaks in order to use them with bedtools-merge
    run: ../../tools/linux-sort.cwl
    in:
      unsorted_file: call_peaks/narrow_peak_file                      # Called peaks as not sorted narrowpeak file
      key:
        default: ["1,1","2,2n","3,3n"]
      output_filename:                                                # Sets the output file name to include narrowpeak file root name and updated suffix
        source: call_peaks/narrow_peak_file
        valueFrom: $(get_root(self.basename)+"_sorted.narrowPeak")
    out:
    - sorted_file                                                     # Coordinate sorted peaks as narrowpeak file

  merge_peaks:
    doc: |
      Merges sorted peaks
    run: ../../tools/bedtools-merge.cwl
    in:
      bed_file: sort_peaks/sorted_file                                # Coordinate sorted peaks as narrowpeak file
      output_filename:                                                # Sets the output file name to include narrowpeak file root name and updated suffix
        source: call_peaks/narrow_peak_file
        valueFrom: $(get_root(self.basename)+"_merged.narrowPeak")
    out:
    - merged_bed_file                                                 # Merged unsorted peaks as narrowpeak file

  sort_merged_peaks:
    doc: |
      Sorts mrged peaks
    run: ../../tools/linux-sort.cwl
    in:
      unsorted_file: merge_peaks/merged_bed_file                      # Merged unsorted peaks as narrowpeak file
      key:
        default: ["1,1","2,2n","3,3n"]
      output_filename:                                                # Sets the output file name to include narrowpeak file root name and updated suffix
        source: merge_peaks/merged_bed_file
        valueFrom: $(get_root(self.basename)+"_sorted.narrowPeak")
    out:
    - sorted_file                                                     # Merged coordinate sorted peaks as narrowpeak file

  count_tags:
    doc: |
      Calculates tag counts for Tn5 binding sites within
      merged sorted peaks
    run: ../../tools/bedtools-intersect.cwl
    in:
      file_a: sort_merged_peaks/sorted_file                           # Merged coordinate sorted peaks as narrowpeak file
      file_b: sort_bed/sorted_file                                    # Coordinate sorted Tn5 binding sites without blacklisted regions
      count:                                                          # For each entry in file_a report the number of hits in file_b or 0
        default: true
      output_filename:                                                # Sets the output file name to include narrowpeak file root name and updated suffix
        source: sort_merged_peaks/sorted_file
        valueFrom: $(get_root(self.basename)+"_tag_counts.bed")
    out:
    - intersected_file

  get_sequences:
    doc: |
      Returns sequences from the merged sorted peaks
    run: ../../tools/bedtools-getfasta.cwl
    in:
      genome_fasta_file: genome_fasta_file
      intervals_file: sort_merged_peaks/sorted_file                   # Merged coordinate sorted peaks as narrowpeak file
    out:
    - sequences_file                                                  # Sequences as fasta file