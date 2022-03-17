cwlVersion: v1.1
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement


inputs:

  indices_folder:
    type: Directory
    doc: "Path to the STAR indices folder"

  fastq_file_r1:
    type: File
    doc: "FASTQ file (read 1). Optionally compressed"

  fastq_file_r2:
    type: File
    doc: "FASTQ file (read 2). Optionally compressed"

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
    default: 1
    doc: "Maximum number of loci the read is allowed to map to"

  max_multimap_anchor:
    type: int?
    default: 50
    doc: "Maximum number of loci anchors are allowed to map to"

  threads:
    type: int?
    default: 1
    doc: "Number of threads for those steps that support multithreading"


outputs:

  qc_report_file_fastq_r1:
    type: File
    outputSource: qc_fastq_r1/html_file
    doc: "FastqQC report for FASTQ file (read 1)"

  qc_report_file_fastq_r2:
    type: File
    outputSource: qc_fastq_r2/html_file
    doc: "FastqQC report for FASTQ file (read 2)"

  trim_adapters_report_r1:
    type: File
    outputSource: trim_adapters/report_file

  trim_adapters_report_r2:
    type: File
    outputSource: trim_adapters/report_file_pair
    doc: "Adapter trimmimg report for FASTQ file (read 2)"

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

  genome_coverage_file:
    type: File
    doc: "Genome coverage file in bigWig format"
    outputSource: calculate_genome_coverage/bigwig_file


steps:

  extract_fastq_r1:
    run: ../../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_r1
    out:
    - fastq_file

  qc_fastq_r1:
    run: ../../tools/fastqc.cwl
    in:
      reads_file: extract_fastq_r1/fastq_file
      threads: threads
    out:
    - html_file

  extract_fastq_r2:
    run: ../../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_r2
    out:
    - fastq_file

  qc_fastq_r2:
    run: ../../tools/fastqc.cwl
    in:
      reads_file: extract_fastq_r2/fastq_file
      threads: threads
    out:
    - html_file

  trim_adapters:
    run: ../../tools/trimgalore.cwl
    in:
      input_file: extract_fastq_r1/fastq_file
      input_file_pair: extract_fastq_r2/fastq_file
      dont_gzip:
        default: true
      length:
        default: 30
      trim1:
        default: true
      paired:
        default: true
    out:
    - trimmed_file
    - trimmed_file_pair
    - report_file
    - report_file_pair

  bypass_trim_adapters:
    run: ../../tools/bypass-trimgalore-pe.cwl
    in:
      original_fastq_file_1: extract_fastq_r1/fastq_file
      trimmed_fastq_file_1: trim_adapters/trimmed_file
      trimming_report_file_1: trim_adapters/report_file
      original_fastq_file_2: extract_fastq_r2/fastq_file
      trimmed_fastq_file_2: trim_adapters/trimmed_file_pair
      trimming_report_file_2: trim_adapters/report_file_pair
      min_reads_count:
        default: 100                                                # any small number should be good
    out:
    - selected_fastq_file_1
    - selected_fastq_file_2

  align_reads:
    run: ../../tools/star-alignreads.cwl
    in:
      readFilesIn:
      - bypass_trim_adapters/selected_fastq_file_1
      - bypass_trim_adapters/selected_fastq_file_2
      genomeDir: indices_folder
      outFilterMultimapNmax: max_multimap
      winAnchorMultimapNmax: max_multimap_anchor
      outFilterMismatchNmax:
        default: 5
      alignSJDBoverhangMin:
        default: 1
      seedSearchStartLmax:
        default: 15
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

  calculate_genome_coverage:
    run: ../../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: sort_alignments_file/bam_bai_pair
      chrom_length_file: select_chrom_length_file/selected_file
      mapped_reads_number:
        source: align_reads/uniquely_mapped_reads_number
        valueFrom: $(self*2)
    out:
    - bigwig_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "RNA Splicing Analysis Paired-End"
label: "RNA Splicing Analysis Paired-End"
s:alternateName: "RNA Splicing Analysis Paired-End"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/devel/rna-splicing-pe.cwl
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
  RNA Splicing Analysis Paired-End
