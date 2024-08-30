# graphical Brain Association Tool (gBAT)

<img align="center" width="1800" height="350" src="https://i.imgur.com/zzfvc5i.jpeg">

The Graphical Brain Association Tool (gBAT) is a MATLAB-based toolbox designed to map brain-behavior associations, with a particular focus on regional mediation analysis. Its purpose is to provide an easy-to-use graphical interface that simplifies exploring the trove of data generated in neuroimaging studies through powerful modeling techniques that can explain relationships between variables through their influence on brain structure or function.

While gBAT emphasizes mediation analysis, it is a highly versatile tool that can be applied to a wide range of brain-behavior analyses. To complement mediation results, the toolbox also supports independent analyses such as correlation, partial correlation, distance correlation, and partial distance correlation. These techniques can help determine whether mediation analysis is appropriate by testing for significant nonlinear relationships between variables. Additionally, these methods can be used on their own, allowing users to explore associations more broadly, all within gBAT’s user-friendly interface.

Core features
•	Graphical interface for mediation and other brain-behavior analyses
•	Automatically aligns brain and behavioral data across diverse naming conventions
•	Automatically bins voxelwise data for atlas-based analysis
•	Automatically plots volumetric brain images (with regional outlines for atlas-based analyses) to visualize results
•	Provides pipeline for automatically projecting results onto fsaverage surface
•	Only implementation of partial distance correlation in matlab, providing identical results to Python’s dcor 0.6 using a different, recursive, formula: https://pypi.org/project/dcor/

To get started, see the pdf manual in this repository. 

gBAT works for voxelwise analyses, but note that we have mostly used it for regional analyses, so if you encounter any problems, please open an issue or reach out to alex.teghipco@sc.edu
