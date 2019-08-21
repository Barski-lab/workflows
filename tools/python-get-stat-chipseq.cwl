cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_filename = function() {
        if (inputs.output_filename == null){
          var root = inputs.bowtie_log.basename.split('.').slice(0,-1).join('.');
          var ext = ".stat";
          return (root == "")?inputs.bowtie_log.basename+ext:root+ext;
        } else {
          return inputs.output_filename;
        }
    };
  - var get_formatted_output_filename = function() {
        if (inputs.formatted_output_filename == null){
          var root = inputs.bowtie_log.basename.split('.').slice(0,-1).join('.');
          var ext = "_formatted.tsv";
          return (root == "")?inputs.bowtie_log.basename+ext:root+ext;
        } else {
          return inputs.formatted_output_filename;
        }
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2


inputs:

  script:
    type: string?
    default: |
      #!/usr/bin/env python
      import sys, re
      TOTAL, ALIGNED, SUPPRESSED, USED = 100, 80, 0, 0
      with open(sys.argv[1], 'r') as bowtie_log:
        for line in bowtie_log:
          if 'processed:' in line:
            TOTAL = int(line.split('processed:')[1])
          if 'alignment:' in line:
            ALIGNED = int(line.split('alignment:')[1].split()[0])
          if 'due to -m:' in line:
            SUPPRESSED = int(line.split('due to -m:')[1].split()[0])
      USED = ALIGNED
      with open(sys.argv[2], 'r') as rmdup_log:
        for line in rmdup_log:
          if '/' in line and 'Skip' not in line:
            splt = line.split('/')
            USED = int((splt[1].split('='))[0].strip()) - int((splt[0].split(']'))[1].strip())
      print TOTAL, ALIGNED, SUPPRESSED, USED
      print >> sys.stderr, "Total reads number\tUniquely mapped reads number\tMulti-mapped reads number\tReads number after removing duplicates"
      print >> sys.stderr, str(TOTAL) + "\t" + str(ALIGNED) + "\t" + str(SUPPRESSED) + "\t" + str(USED)

    inputBinding:
      position: 5
    doc: |
      Python script to get TOTAL, ALIGNED, SUPPRESSED, USED values from log files

  bowtie_log:
    type: File
    inputBinding:
      position: 6
    doc: |
      Log file from Bowtie

  rmdup_log:
    type: File
    inputBinding:
      position: 7
    doc: |
      Log file from samtools rmdup

  output_filename:
    type:
    - "null"
    - string
    doc: |
      Name for generated output file

  formatted_output_filename:
    type:
    - "null"
    - string
    doc: |
      Name for generated formatted output file

outputs:

  output_file:
    type: File
    outputBinding:
      glob: $(get_output_filename())

  formatted_output_file:
    type: File
    outputBinding:
      glob: $(get_formatted_output_filename())

  total_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[0]))

  mapped_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[1]))

  supressed_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[2]))

  used_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[3]))

baseCommand: [python, '-c']
arguments:
  - valueFrom: $(" > " + get_output_filename() + " 2> " + get_formatted_output_filename())
    position: 100000
    shellQuote: false

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "python-get-stat-chipseq"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/python-get-stat-chipseq.cwl
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
  Tool processes and combines log files generated by Bowtie aligner and samtools rmdup.

  `get_output_filename` function returns output filename equal to `output_filename` (if this input is provided) or
  generated on the base of bowtie log basename with `.stat` extension.

  `get_formatted_output_filename` function returns output filename equal to `formatted_output_filename` (if input is provided) or
  generated on the base of STAR log basename with `_formatted.tsv` extension.

s:about: |
  Runs python code from the script input
