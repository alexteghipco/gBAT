# basicVBM
These scripts may help you perform some basic voxel based morphometry with correlation, partial correlation, distance correlation, partial distance correlation or mediation analysis. They also support correction for cluster size with permutations and ROI-based analyses. See demo_start_here.m to get started (new demo, will add more details in a forthcoming update). Note, you will have to plug in and organize your own behavioral and imaging data for the demo for now but we do use a publicly available dataset and at some point will share a link to download the demo data.

Dependencies
Please download the following packages/data and make sure they are in your MATLAB path: 
  
  1. Mediation toolbox from Tor Wager's lab: https://github.com/canlab/MediationToolbox
      * we will rely on just  mediation.m and its associated functions
  2. Auxillary functions for the mediation toolbox, which are part of the larger lab repo: https://github.com/canlab/CanlabCore
  3. Check to make sure that you have the bioinformatics and image processing toolbox from matlab installed
     * type ver into the command window, hit enter, and search for the toolbox
     * we will use these for FDR correction and finding connected components (i.e., clusters) in images
  4. If you plan to use distance correlation or partial distance correlation for VBM, download this package and make sure it is in your path: https://github.com/alexteghipco/partialDistanceCorrelation
  5. The toolbox comes with a function for loading in brain data into matlab but you will need either SPM, matlab's image processing toolbox or freesurfer's matlab toolbox in your path




