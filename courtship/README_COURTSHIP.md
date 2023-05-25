# Courtship Song Analysis README
MATLAB code used to analyze and plot data from courtship song experiments. For details of experimental protocol, see Methods section of manuscript. 

The bulk of the analysis comes from code written by the Stern and Murthy Labs:
  * **FlySongSegmenter**: automatic segmentation of *Drosophila* courtship song (Paper: [Arthur et al., 2013](https://doi.org/10.1186/1741-7007-11-11) | GitHub: [FlySongSegmenter](https://github.com/FlyCourtship/FlySongSegmenter.git))
  * **BatchSongAnalysis**: bulk statistical analysis of segmented *Drosophila* courtship song data (Stern Lab | GitHub: [BatchSongAnalysis](https://github.com/dstern/BatchSongAnalysis.git))
  * **pulseTypeClassifier**: classification of song pulses into fast and slow sub-types (Murthy Lab | Paper: [Clements et al., 2014](https://doi.org/10.1038/nature13131) | GitHub: [pulseTypeClassifier](https://github.com/postpop/pulseTypeClassifier.git)). NB: While the linked pulseTypeClassifier repository is a self-contained version of the pulse classifier, the code in *this* repository was taken from the [MurthyLab_FlySongSegmenter](https://github.com/murthylab/MurthyLab_FlySongSegmenter.git) to include the additional 'tools' folder. 

## Directory structure
This directory is organized as follows:
  * *analysis/*: code which, alongside the directories in the *external/* folder, comprises the scripts and functions used in the data analysis pipeline.
  * *data/*: data structures necessary to reproduce courtship song figures. Also contains *example.zip* which, when extracted, allows testing of the analysis pipeline that is not covered in other GitHub repositories--the data in *example.zip* is what one would get after having run `BatchSongAnalysis.m`
  * *external/*: external code sources used in the analysis pipeline (see above).
  * *figs/*: code to produce courtship song figures. The two scripts `plotSongStats.m` and `plotPulseFig.m` produce the two components of the courtship figure in the manuscript. Also contains:
	  * *output/*: a folder to save figure output
	  * *subtightplot/*: the MATLAB subtightplot tool used for arranging plot axes (Felipe G. Nievinski | MATLAB File Exchange: [subtightplot](https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot))

## Analysis pipeline
For analysis of proportion of time spent singing, as well as proportion of song time devoted to pulse vs sine song:
  1. Data collected from courtship song experiments, in the form of .daq files, is segmented using `external/FlySongSegmenter/FlySongSegmenter.m`
  2. Batch analysis of segmented songs is performed using `external/BatchSongAnalysis/BatchFlySongAnalysis.m`
  3. File paths and genotype information is collected by `analysis/getDataDirs.m`, which creates a dataPathStruct containing the aforementioned info (NB: this requires the user to enter some path information; may require some fiddling to get to work with user's preferred directory structure for data storage)
  4. Time spent singing, sine index, and pulse index are aggregated with `analysis/getSongStats.m`

For pulse classification (assumes song has been segmented/analyzed above already):
  1. Classification of pulses performed across flies by `analysis/batchPulseClassifier.m` (this script calls the classifier from Murthy Lab pulseTypeClassifier noted above)
  2. Data from classification analysis across genotypes is aggregated by `analysis/makePulseSummaryStruct.m`
