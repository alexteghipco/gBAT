# basicVBM

<img align="center" width="800" height="800" src="https://i.imgur.com/iS3jdzb.png">

This is a tool with a GUI that might be able to help you perform some basic voxel based morphometry with correlation, partial correlation, distance correlation, partial distance correlation or mediation analysis (or to use these methods to relate any brain data to some behavioral measures). Technically, you could use this to look for relationships in any brain data but it has not used this way yet.

To start, download this repository, open matlab, navigate to the your downloaded folder, and type vbmGui into the command window. This will summon the GUI. Please hover over buttons to see detailed information about what they do and how to use them.

Some notes about usage. The GUI disables some buttons if you don't have the functions required to make them work. Depending on the kind of analysis you want to do, you'll need some of these other toolboxes in your path: 
  1. For mediation analysis, the Mediation toolbox from Tor Wager's lab: https://github.com/canlab/MediationToolbox (we use their excellent mediation.m function to perform the mediation analysis and someetimes do our own bootstrapping for significance)
  2. You will need their bigger lab repo if you want to generate path diagrams and partial regression plots (to do so, add plot to the additional input argument section of the GUI):
     https://github.com/canlab/CanlabCore

You may need some other MATLAB toolboxes depending on your use-case. For example, statistics and machine learning for correlation and partial correlation, Bioinformatics for FDR correction (always the original method), and imaging toolbox for loading in and writing nifti files to analyze. Note, we repackage our partial distance correlation repo with this. This package has a unique implementation for computing partial distance correlation (see readme in that folder) but we use the very nice bias corrected distance correlation function made available here: https://www.mathworks.com/matlabcentral/fileexchange/58445-bias-corrected-distance-correlation within that package. Some colormap functions are included with this repo as well in order to make brain montages. These have also not been written by me. Please consider citing the wonderful open source work of all of these other authors if you use any of their code. Relevant information can be found within the functions and this readme will be updated with relevant citation information to streamline this process.

vbmGui generates lots of files for you to look at that will be self-explanatory (but more info forthcoming). If you want to apply a new p-value threshold, you can use the outputs in your folder to generate new maps (provided you used FDR correction or no correction). To do this, see the supplemental files in ./postpro (note these have not been debugged as thoroughly). This will also show you how to generate volume space brain images like vbmGui automatically does as part of its process. Finally, you can use atlasVol2Surf to copy a volume space result for ROIs within some atlas onto an existing surface space version of that atlas, converting your ROI result into surface space without having to individually project ROIs, etc.

Some warnings:
  1. The GUI has been debugged for ROI-based analyses. Please see ./vbmHelpers for how your atlas information should be set up (see
     rjhu.nii.gz and jhu.txt)
  2. vbmGui tries to align your behavioral and brain data automatically but naming conventions vary between labs and projects. Please
     ensure correct mapping (there will be some warnings that pop up). The attempt may fail as it just looks for the first column of your
     behavioral data and tries to find nifti files that contain that same string. In ties, it will go with the shortest file name.
  4. Voxelwise analyses still need to be debugged in the GUI but you can use the old code from the demo_start_here.m to setup your own
     analysis without a GUI. This should work, though there have been some changes to the structure of vbm.m that might break things (update
     forthcoming, use the older version of this repo if you have it for now).
  6. If you ask vbmGui to generate brain figures it will take a bit of time. We use matlab's contour function to outline and plot each ROI
     individually at the moment so it is a laborious process. Brain images are NOT generated using matlab's montager function as it lacks
     some critical functionality.
  8. Dependening on your analysis you may need a lot of RAM (for example, distance correlation, where we need to get euclidean distances
     between all samples for many variables, but also mediation or if bootstrapping for any analysis is enabled).
  9. An option is available to paralelize models being generated in vbm.m but this may backfire and make things slower if your data is high
      dimensional (lots or ROIs or voxels)
