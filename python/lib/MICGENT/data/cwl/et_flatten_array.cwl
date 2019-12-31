cwlVersion: v1.0
class: ExpressionTool
doc: Given x[a[],b[],c[],...], return a[]+b[]+c[]+... where `+` means concatenation. Use to flatten output from a scatter.
requirements:
  - class: InlineJavascriptRequirement
inputs:
  inp: Any[]
outputs:
  flat: File[]
expression: |
  ${
    return {"flat": [].concat.apply([], inputs.inp)};
  }