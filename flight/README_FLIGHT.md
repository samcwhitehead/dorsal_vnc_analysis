# Flight Analysis README
MATLAB code used to analyze and plot data from tethered flight experiments collected using [Kinefly](https://github.com/ssafarik/Kinefly.git). For detailed methods of experimental protocol, see Methods section of manuscript. 

## Directory structure
This directory is organized as follows:
  * *analysis/*: code which comprises the scripts and functions used in the data analysis pipeline.
  * *data/*: data structures necessary to reproduce tethered flight figures. 
  * *figs/*: code to produce tethered flight analysis figures. `make_example_fly_plot.m` and `make_kf_summary_plot.m` create the figure panels found in the manuscript. Folder also contains the subfolder *output/*, where created figures can be saved.

## Analysis pipeline
1. .abf files are obtained directly from tethered flight experiments, and contain the stored kinematics tracked by Kinefly
2. Raw data (.abf files) are read and converted to .mat files using `analysis/run_kinefly_analysis.m`
3. Data from multiple genotypes/experimental conditions are aggregated using `analysis/combine_kinefly_data.m` which produces two types of .mat file: i) a time_series_struct which contains the full temporal resolution data from the .abf file for all flies from a single data folder and ii) a grand_mean_struct which contains sub-sampled data from multiple folders as well as a grand mean and bootstrapped confidence intervals for each of the kinematics measurements (right and left wing amplitude, wingbeat frequency, etc.). NB: the code assumes that data is organized as follows: dataRoot/genotype/condition/fly/filename.abf (e.g. ../SS41039/closed_loop_50x/20190210T165714_SS41039/blah.abf). If your data is not organized this way you may need to tweak the code as written.
4. Results are plotted using the code in *figs/*
