#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.3-2
# To generate again: $ ariba --generate_cwl_tool
# Help: $ ariba --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['ariba', 'micplot']

doc: |
  Make violin/dot plots using MIC data

inputs:
  
  prepareref_dir:
    type: string
  
    doc: Name of output directory when "ariba prepareref" was run
    inputBinding:
      position: 1

  antibiotic:
    type: string
  
    doc: Antibiotic name. Must exactly match a column from the MIC file
    inputBinding:
      position: 2

  mic_file:
    type: string
  
    doc: File containing MIC data for each sample and one or more antibiotics
    inputBinding:
      position: 3

  summary_file:
    type: string
  
    doc: File made by running "ariba summary"
    inputBinding:
      position: 4

  outprefix:
    type: string
  
    doc: Prefix of output files
    inputBinding:
      position: 5

  out_format:
    type: ["null", string]
    default: pdf
    doc: Output format of image file. Use anything that matplotlib can save to, eg pdf or png [%(default)s]
    inputBinding:
      prefix: --out_format 

  main_title:
    type: ["null", string]
    doc: Main title of plot. Default is to use the antibiotic name
    inputBinding:
      prefix: --main_title 

  plot_height:
    type: ["null", float]
    default: 7
    doc: Height of plot in inches [%(default)s]
    inputBinding:
      prefix: --plot_height 

  plot_width:
    type: ["null", float]
    default: 7
    doc: Width of plot in inches [%(default)s]
    inputBinding:
      prefix: --plot_width 

  use_hets:
    type:
    - "null"
    - type: enum
      symbols: ['yes', 'no', 'exclude']
    default: yes
    doc: How to deal with HET snps. Choose from yes,no,exclude. yes - count a het SNP as present. no - do not count a het SNP as present. exclude - completely remove any sample with any het SNP [%(default)s]
    inputBinding:
      prefix: --use_hets 

  interrupted:
    type: ["null", boolean]
    default: False
    doc: Include interrupted genes (as in the assembled column of the ariba summary files)
    inputBinding:
      prefix: --interrupted 

  min_samples:
    type: ["null", int]
    default: 1
    doc: Minimum number of samples in each column required to include in plot [%(default)s]
    inputBinding:
      prefix: --min_samples 

  no_combinations:
    type: ["null", boolean]
    default: False
    doc: Do not show combinations of variants. Instead separate out into one box/violin plot per variant.
    inputBinding:
      prefix: --no_combinations 

  panel_heights:
    type: ["null", string]
    default: 9,2
    doc: Two integers that determine relative height of top and bottom plots. eg 5,1 means ratio of 5 -1 between top and bottom panel heights [%(default)s]
    inputBinding:
      prefix: --panel_heights 

  panel_widths:
    type: ["null", string]
    default: "5,1"
    doc: Two integers that determine relative width of plots and space used by counts legend. eg 5,1 means ratio of 5 -1 between top and bottom panel widths. Only applies when plotting points and --point_size 0 [%(default)s]
    inputBinding:
      prefix: --panel_widths 

  count_legend_x:
    type: ["null", float]
    default: -2
    doc: Control x position of counts legend when plotting points and --point_size 0 [%(default)s]
    inputBinding:
      prefix: --count_legend_x 

  p_cutoff:
    type: ["null", float]
    default: 0.05
    doc: p-value cutoff for Mann-Whitney tests [%(default)s]
    inputBinding:
      prefix: --p_cutoff 

  xkcd:
    type: ["null", boolean]
    default: False
    doc: Best used with xkcd font installed ;)
    inputBinding:
      prefix: --xkcd 

  colourmap:
    type: ["null", string]
    default: "Accent"
    doc: Colours to use. See http -//matplotlib.org/users/colormaps.html [%(default)s]
    inputBinding:
      prefix: --colourmap 

  number_of_colours:
    type: ["null", int]
    default: 0
    doc: Number of colours in plot. 0 -same number as columns in the plot. 1 -all black. >1 - take the first N colours from the colourmap specified by --colourmap and cycle them [%(default)s]
    inputBinding:
      prefix: --number_of_colours 

  colour_skip:
    type: ["null", string]
    doc: If using a continuous colourmap, --colour_skip a,b (where 0 <= a < b <= 1) will skip the range between a and b. Useful for excluding near-white colours
    inputBinding:
      prefix: --colour_skip 

  plot_types:
    type: ["null", string]
    default: "violin,point"
    doc: Types of plots to make, separated by commas. Choose from violin,point [%(default)s]
    inputBinding:
      prefix: --plot_types 

  hlines:
    type: ["null", string]
    default: 
    doc: Comma-separated list of positions at which to draw horizontal lines. Default is to draw no lines.
    inputBinding:
      prefix: --hlines 

  jitter_width:
    type: ["null", float]
    default: 0.1
    doc: Jitter width option when plotting points [%(default)s]
    inputBinding:
      prefix: --jitter_width 

  log_y:
    type: ["null", float]
    default: 2
    doc: Base of log applied to y values. Set to zero to not log [%(default)s]
    inputBinding:
      prefix: --log_y 

  point_size:
    type: ["null", float]
    default: 4
    doc: Size of points when --plot_types includes point. If zero, will group points and size them proportional to the group size [%(default)s]
    inputBinding:
      prefix: --point_size 

  point_scale:
    type: ["null", float]
    default: 1
    doc: Scale point sizes when --point_size 0. All point sizes are multiplied by this number. Useful if you have large data set [%(default)s]
    inputBinding:
      prefix: --point_scale 

  violin_width:
    type: ["null", float]
    default: 0.75
    doc: Width of violins [%(default)s]
    inputBinding:
      prefix: --violin_width 

  dot_size:
    type: ["null", float]
    default: 100
    doc: Size of dots in lower part of plot [%(default)s]
    inputBinding:
      prefix: --dot_size 

  dot_outline:
    type: ["null", boolean]
    default: False
    doc: Black outline around all dots (whether coloured or not) in lower part of plots
    inputBinding:
      prefix: --dot_outline 

  dot_y_text_size:
    type: ["null", int]
    default: 7
    doc: Text size of labels [%(default)s]
    inputBinding:
      prefix: --dot_y_text_size 


outputs:
    []
