cwlVersion: v1.1
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_root = function(basename) {
          return basename.split('.').slice(0,1).join('.');
      };


inputs:

  indices_folder:
    type: Directory
    doc: "Path to the STAR indices folder"

  fastq_file:
    type: File
    doc: "FASTQ file. Optionally compressed"

  clip_3p_nbases:
    type: int?
    default: 0
    doc: "Number of bases to clip from the 3p end"

  clip_5p_nbases:
    type: int?
    default: 0
    doc: "Number of bases to clip from the 5p end"

  max_multimap:
    type: int?
    default: 10
    doc: "Maximum number of loci the read is allowed to map to"

  max_multimap_anchor:
    type: int?
    default: 50
    doc: "Maximum number of loci anchors are allowed to map to"

  out_filter_mismatch_nmax:
    type: int?
    default: 10
    doc: "Alignment will be output only if it has no more mismatches than this value"

  align_sjdb_overhang_min:
    type: int?
    default: 3
    doc: "Minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments"

  seed_search_start_lmax:
    type: int?
    default: 50 
    doc: "Defines the search start point through the read - the read is split into pieces no longer than this value"

  threads:
    type: int?
    default: 1
    doc: "Number of threads for those steps that support multithreading"


outputs:

  qc_report_file_fastq:
    type: File
    outputSource: qc_fastq/html_file
    doc: "FastqQC report for FASTQ file"

  trim_adapters_report:
    type: File
    outputSource: trim_adapters/report_file

  align_reads_final_log:
    type: File?
    outputSource: align_reads/log_final
    doc: "STAR Log.final.out"

  align_reads_out_log:
    type: File?
    outputSource: align_reads/log_out
    doc: "STAR Log.out"

  align_reads_progress_log:
    type: File?
    outputSource: align_reads/log_progress
    doc: "STAR Log.progress.out"

  align_reads_stdout_log:
    type: File?
    outputSource: align_reads/log_std
    doc: "STAR Log.std.out"

  align_reads_sj_log:
    type: File?
    outputSource: align_reads/log_sj
    doc: "STAR SJ.out.tab"

  alignments_file:
    type: File
    outputSource: sort_alignments_file/bam_bai_pair
    doc: "Coordinate sorted BAM file and BAI index file"

  qc_report_file_alignments:
    type: File
    outputSource: qc_alignments_file/log_file
    doc: "BAM file statistics report"

  genome_coverage_not_stranded_file:
    type: File
    doc: "Not stranded genome coverage file in bigWig format"
    outputSource: calculate_genome_coverage_not_stranded/bigwig_file

  genome_coverage_positive_file:
    type: File
    doc: "Genome coverage file for positive strand in bigWig format"
    outputSource: calculate_genome_coverage_positive/bigwig_file

  genome_coverage_negative_file:
    type: File
    doc: "Genome coverage file for negative strand in bigWig format"
    outputSource: calculate_genome_coverage_negative/bigwig_file


steps:

  extract_fastq:
    run: ../../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
      output_prefix:
        default: "read_1"
    out:
    - fastq_file

  qc_fastq:
    run: ../../tools/fastqc.cwl
    in:
      reads_file: extract_fastq/fastq_file
      threads: threads
    out:
    - html_file

  trim_adapters:
    run: ../../tools/trimgalore.cwl
    in:
      input_file: extract_fastq/fastq_file
      length:
        default: 30
    out:
    - trimmed_file
    - report_file

  align_reads:
    run: ../../tools/star-alignreads.cwl
    in:
      readFilesIn: trim_adapters/trimmed_file
      genomeDir: indices_folder
      outFilterMultimapNmax: max_multimap
      winAnchorMultimapNmax: max_multimap_anchor
      outFilterMismatchNmax: out_filter_mismatch_nmax
      alignSJDBoverhangMin: align_sjdb_overhang_min
      seedSearchStartLmax: seed_search_start_lmax
      outSAMattributes:
        default: "All"                              # for SplAdder we need NM tag
      twopassMode:
        default: "Basic"
      clip3pNbases: clip_3p_nbases
      clip5pNbases: clip_5p_nbases
      threads: threads
    out:
    - aligned_file
    - log_final
    - uniquely_mapped_reads_number
    - log_out
    - log_progress
    - log_std
    - log_sj

  sort_alignments_file:
    run: ../../tools/samtools-sort-index.cwl
    in:
      sort_input: align_reads/aligned_file
      threads: threads
    out:
    - bam_bai_pair

  qc_alignments_file:
    run: ../../tools/samtools-stats.cwl
    in:
      bambai_pair: sort_alignments_file/bam_bai_pair
    out:
    - log_file

  select_chrom_length_file:
    run: ../../tools/get-file-by-name.cwl
    in:
      input_files: indices_folder
      basename_regex:
        default: "chrNameLength.txt"
    out:
    - selected_file

  calculate_genome_coverage_not_stranded:
    run: ../../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: sort_alignments_file/bam_bai_pair
      chrom_length_file: select_chrom_length_file/selected_file
      mapped_reads_number:
        source: align_reads/uniquely_mapped_reads_number
        valueFrom: $(self*2)
      bigwig_filename:
        source: sort_alignments_file/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_not_stranded.bigWig")
    out:
    - bigwig_file

  calculate_genome_coverage_positive:
    run: ../../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: sort_alignments_file/bam_bai_pair
      chrom_length_file: select_chrom_length_file/selected_file
      mapped_reads_number:
        source: align_reads/uniquely_mapped_reads_number
        valueFrom: $(self*2)
      bigwig_filename:
        source: sort_alignments_file/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_positive.bigWig")
      strand:
        default: "+"
      dutp:
        default: true
    out:
    - bigwig_file

  calculate_genome_coverage_negative:
    run: ../../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: sort_alignments_file/bam_bai_pair
      chrom_length_file: select_chrom_length_file/selected_file
      mapped_reads_number:
        source: align_reads/uniquely_mapped_reads_number
        valueFrom: $(self*2)
      bigwig_filename:
        source: sort_alignments_file/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_negative.bigWig")
      strand:
        default: "-"
      dutp:
        default: true
    out:
    - bigwig_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "RNA Splicing Analysis Single Read"
label: "RNA Splicing Analysis Single Read"
s:alternateName: "RNA Splicing Analysis Single Read"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/devel/rna-splicing-se.cwl
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
  RNA Splicing Analysis Single Read
