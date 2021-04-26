cwlVersion: v1.0
class: CommandLineTool

requirements: []

hints:
  - class: DockerRequirement
    dockerPull: haysb1991/altanalyze-test:version9

inputs:
  move_script:
    type: string?
    default: |
      #!/bin/bash
      echo "$@"
      mkdir /opt/altanalyze/userdata/
      cp -r "$0" /opt/altanalyze/AltDatabase
      cp "$1" /opt/altanalyze/userdata/
      cp "$2" /opt/altanalyze/userdata/
      cp "$3" /opt/altanalyze/userdata/
      cp "$4" /opt/altanalyze/userdata/
      ls /opt/altanalyze/userdata
      python /opt/altanalyze/AltAnalyze.py --cellHarmony yes --species Mm --platform RNASeq --fold 1.2 --adjp yes --pval 0.05 --correlationCutoff 0.7 --referenceType centroid --reference "$1" --labels "$2" --referenceFull "$3" --input "$4"
      ls /opt/altanalyze/userdata
      cp -r /opt/altanalyze/userdata .

    inputBinding:
      position: 1
    doc: |
      Bash function to redirect to complete the return of EnsMart72 as output.

  data_in:
    type: Directory
    inputBinding:
      position: 5

  reference:
    type: File
    inputBinding:
      position: 6

  labels:
    type: File
    inputBinding:
      position: 7

  referenceFull:
    type: File
    inputBinding:
      position: 8

  input:
    type: File
    inputBinding:
      position: 9

outputs:
  stdout_log:
    type: stdout

  stderr_log:
    type: stderr

  ICGSNMF:
    type: Directory
    outputBinding: 
      glob: "userdata"

baseCommand: [bash, '-c']

stdout: chtest_stdout.log
stderr: chtest_stderr.log

$namespaces:
  s: http://schema.org/

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

doc: |
  altanalyze is being used in a single cell pipeline