cwlVersion: v1.0
class: CommandLineTool

requirements: []

hints:
  - class: DockerRequirement
    dockerPull: haysb1991/python3.8_cellbrowser

inputs:
  move_script:
    type: string?
    default: |
      #!/bin/bash
      echo "$@"
      mkdir /outputhtml/
      cp "$0" exprMatrix.tsv
      cp "$1" meta.tsv
      cp "$2" tsne.coords.tsv
      cp "$3" .
      cp "$4" .
      cp "$5" .
      cbBuild -o outputhtml
      ls outputhtml
      cp -r outputhtml .

    inputBinding:
      position: 1
    doc: |
      Bash function to generate needed html and scripts for cellbrowser.

  input1:
    type: File
    inputBinding:
      position: 2

  input2:
    type: File
    inputBinding:
      position: 3

  input3:
    type: File
    inputBinding:
      position: 4

  input4:
    type: File
    inputBinding:
      position: 5

  input5:
    type: File
    inputBinding:
      position: 6

  input6:
    type: File
    inputBinding:
      position: 7

outputs:
  stdout_log:
    type: stdout

  stderr_log:
    type: stderr

  htmlout:
    type: Directory
    outputBinding:
      glob: "outputhtml"    

baseCommand: [bash, '-c']

stdout: setupcb_stdout.log
stderr: setupcb_stderr.log

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
  cellbrowser is an interactive visualization tool