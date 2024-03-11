# basicVBM
This is a tool with a GUI that might be able to help you perform some basic voxel based morphometry with correlation, partial correlation, distance correlation, partial distance correlation or mediation analysis (or to use these methods to relate any brain data to some behavioral measures). 

To start, download this repository, open matlab, navigate to the your downloaded folder, and type vbmGui into the command window. This will summon the GUI. Please hover over buttons to see detailed information about what they do and how to use them.

Some notes about usage. The GUI disables some buttons if you don't have the functions required to make them work. Depending on the kind of analysis you want to do, you'll need some of these other toolboxes in your path: 
  1. For mediation analysis, the Mediation toolbox from Tor Wager's lab: https://github.com/canlab/MediationToolbox (we use their excellent      mediation.m function to perform the mediation analysis and someetimes do our own bootstrapping for significance)
  2. You will need their bigger lab repo if you want to generate path diagrams and partial regression plots:
     https://github.com/canlab/CanlabCore

You may need some other MATLAB toolboxes depending on your use-case. For example, statistics and machine learning for correlation and partial correlation, Bioinformatics for FDR correction (always the original method), and imaging toolbox for loading in and writing nifti files to analyze. 

vbmGui generates lots of files for you to look at that will be self-explanatory. If you want to apply a new p-value threhsold, you can use the outputs in your folder to generate new maps (provided you used FDR correction or no correction). To do this, see the supplemental files in ./postpro. This will also show you how to generate volume space brain images like vbmGui automatically does as part of its process. Finally, you can use atlasVol2Surf to copy a volume space result for ROIs within some atlas onto an existing surface space version of that atlas, converting your ROI result into surface space without having to individually project ROIs, etc.

Some warnings:
  1. The GUI has been debugged for ROI-based analyses. Please see ./vbmHelpers for how your atlas information should be set up (see
     rjhu.nii.gz and jhu.txt)
  2. vbmGui tries to align your behavioral and brain data automatically but naming conventions vary between labs and projects. Please
     ensure correct mapping (there will be some warnings that pop up). The attempt may file as it just looks for the first column of your
     behavioral data and tries to find nifti files that contain that same string. In ties, it will go with the shortest file name.
  4. Voxelwise analyses still need to be debugged in the GUI but you can use the old code from the demo_start_here.m to setup your own
     analysis without a GUI. This should work, though there have been some changes to the structure of vbm.m that might break things (update
     forthcoming).
  6. If you ask vbmGui to generate brain figures it will take a bit of time. We use matlab's contour function to outline and plot each ROI
     individually at the moment so it is a laborious process. Brain images are NOT generated using matlab's montager.m function as it lacks
     some critical functionality.
  8. Dependening on your analysis you may need a lot of RAM (for example, distance correlation, where we need to get euclidean distances
     between all samples for many variables, but also mediation or if bootstrapping for any analysis is enabled).
  9. An option is available to paralelize models being generated in vbm.m but this may backfire and make things slower if your data is high      dimensional (lots or ROIs or voxels)
