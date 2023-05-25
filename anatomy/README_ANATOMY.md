# Anatomy Analysis README
MATLAB code used to analyze the anatomy, input/output profile, and hemilineage identity of neurons targeted by transgenic driver lines in the manuscript and segmented from multi-color flip-out (MCFO) confocal microscopy images. Also contains MATLAB code for generating figures based on these analyses. For details of methodology, see Methods section of manuscript. 

## Directory structure
This directory is organized as follows:
  * *analysis/*: code which comprises the scripts and functions used in the data analysis pipeline. This is subdivided into:
	  * *clustering/*: used to perform hierarchical clustering of cells based on the VNC neuropils in which they have inputs/outputs.
	  * *image_processing/*: used to process and binarize segmented confocal stacks of neurons for use in subsequent analyses
	  * *overlap/*: used to calculate masked volume overlap between neuron pairs. 
  * *data/*: data structures necessary to reproduce anatomy figures. 
  * *figs/*: code to produce anatomical analysis figures. This is subdivided into:
	  * *clustering/*: produces cluster matrix and directed graph plots
	  * *misc/*: contains miscellaneous plot tools
	  * *overlap/*: produces plots for masked volume overlap between neurons
  * *util/*: helper code for data analysis/handling.

