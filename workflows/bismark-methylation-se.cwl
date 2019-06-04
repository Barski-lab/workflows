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

  threads:
    type: int?
    default: 1
    label: "Number of cores to use"
    doc: "Sets the number of parallel instances of Bismark to be run concurrently"


outputs:

  bambai_pair:
    type: File
    label: "BAM + BAI files"
    doc: "Bismark aligned coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index/bam_bai_pair

  bismark_align_log_file:
    type: File
    label: "Bismark aligner log file"
    doc: "Log file generated by Bismark on the alignment step"
    outputSource: bismark_align/log_file

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
    outputSource: refactore_mbias_plot/refactored_mbias_plot

  bedgraph_cov_file:
    type: File
    label: "Methylation statuses in bedGraph format"
    doc: "Methylation statuses in bedGraph format"
    outputSource: bismark_extract_methylation/bedgraph_cov_file

  bismark_cov_file:
    type: File
    label: "Genome-wide cytosine methylation report"
    doc: "Coverage text file summarising cytosine methylation values"
    outputSource: bismark_extract_methylation/bismark_cov_file

  splitting_report_file:
    type: File
    label: "Methylation extraction log"
    doc: "Log file giving summary statistics about methylation extraction"
    outputSource: bismark_extract_methylation/splitting_report_file


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
      threads: threads
    out: [bam_file, log_file]

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
      threads: threads
    out:
      - chg_context_file
      - chh_context_file
      - cpg_context_file
      - mbias_plot
      - bedgraph_cov_file
      - bismark_cov_file
      - splitting_report_file

  refactore_mbias_plot:
    in:
      mbias_plot: bismark_extract_methylation/mbias_plot
    out: [refactored_mbias_plot]
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: refactore.py
            entry: |
                #!/usr/bin/env python
                import sys
                prefix_list, i, prefix, lines = ["CpG", "CHG", "CHH"], 0, "", {}
                with open(sys.argv[1], 'r') as infile:
                    for line in infile:
                        if line.strip():
                            if "context" in line:
                                prefix = ""
                            if "position" in line:
                                prefix = prefix_list[i]
                                i += 1
                                continue
                            if not prefix:
                                continue
                            data = line.split()
                            pos = int(data[0])
                            insert_data = {prefix + " count methylated":   data[1],
                                          prefix + " count unmethylated": data[2],
                                          prefix + " % methylation":      data[3] if (data[1] != "0" and data[2] != "0") else 0,
                                          prefix + " coverage":           data[4] if (data[1] != "0" and data[2] != "0") else 0}
                            if pos in lines:
                                lines[pos].update(insert_data)
                            else:
                                lines[pos] = insert_data
                print("positions\tCHG % methylation\tCHG count methylated\tCHG count unmethylated\tCHG coverage\tCHH % methylation\tCHH count methylated\tCHH count unmethylated\tCHH coverage\tCpG % methylation\tCpG count methylated\tCpG count unmethylated\tCpG coverage")
                for k, v in sorted(lines.items()):
                  print str(k), "\t", "\t".join([str(t) for s,t in sorted(v.items())])
      - class: DockerRequirement
        dockerPull: biowardrobe2/python-pandas:v0.0.1
      inputs:
        mbias_plot:
          type: File
          inputBinding:
            position: 1
      outputs:
        refactored_mbias_plot:
          type: stdout
      stdout: "mbias_plot.tsv"
      baseCommand: ["python", "refactore.py"]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "bismark-methylation"
s:downloadUrl: https://github.com/Barski-lab/workflows/blob/master/workflows/bismark-methylation.cwl
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





