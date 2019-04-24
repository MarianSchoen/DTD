# DTD
Digital Tissue Deconvolution (DTD) reconstructs the cellular composition of a tissue from its bulk expression profile. 
In order to increase deconvolution accuracy, DTD adapts the deconvolution model to the tissue scenario via loss-function learning.  
Training is performed on 'in-silicio' training mixtures, for which the cellular composition are known.  
As input, DTD requires a labelled expression matrix. 
The package includes functions to generate training and test mixtures, train the model, and assess its deconvolution capability via visualizations.  
In Goertler et al. 2018 "Loss-function Learning for Digital Tissue Deconvolution" the theory has been published.

An exemplary analysis can be viewed at https://github.com/MarianSchoen/Exemplary-DTD-analysis

# Install
Install from github, without vignette: 
``` r
  devtools::install_github("MarianSchoen/DTD")
```
I strongly recommend creating the vignette. 
Therefore, install from github with vignette  
(creating vignettes approximately takes ~3 minutes)
``` r
  devtools::install_github(
    "MarianSchoen/DTD", 
    build_opts = c("--no-resave-data", "--no-manual"), 
    build_vignettes=TRUE
    )
  browseVignettes("DTD")
```

# Introduction to DTD Theory
The gene expression profile of a tissue combines the expression profiles of all cells in this tissue. Digital tissue deconvolution (DTD) addresses the following inverse problem: Given the expression profile y of a tissue, what is the cellular composition c of cells X in that tissue? The cellular composition c can be estimated by  
  
<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{150}&space;arg&space;~min_c&space;||y&space;-&space;Xc||_2^2" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;arg&space;~min_c&space;||y&space;-&space;Xc||_2^2" title="arg ~min_c ||y - Xc||_2^2" /></a> 
  


Görtler et al (2019) generalized this formula by introducing a vector g  

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{150}&space;arg&space;~min_c&space;||diag(g)(y&space;-&space;Xc)||_2^2~~~~~~~~~~~~~~~~~~~~(2)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;arg&space;~min_c&space;||diag(g)(y&space;-&space;Xc)||_2^2~~~~~~~~~~~~~~~~~~~~(2)" title="arg ~min_c ||diag(g)(y - Xc)||_2^2~~~~~~~~~~~~~~~~~~~~(2)" /></a> 

Every entry g[i] of g holds the information how important gene i is for the deconvolution process. It can either be selected via prior knowledge, or learned on training data. Training data consists of artificial bulk profiles Y, and the corresponding cellular compositions C. We generate this data with single cell RNASeq profiles.  
The underlying idea of loss-function learning DTD is to gain the vector g by minimizing a loss function L on the training set: 

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{150}&space;L&space;=&space;-\sum_j&space;cor(C_{j,&space;.},&space;\widehat{C_{j,&space;.}}(g))&space;&plus;&space;\lambda&space;||g||_1" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;L&space;=&space;-\sum_j&space;cor(C_{j,&space;.},&space;\widehat{C_{j,&space;.}}(g))&space;&plus;&space;\lambda&space;||g||_1" title="L = -\sum_j cor(C_{j, .}, \widehat{C_{j, .}}(g)) + \lambda ||g||_1" /></a>

Here, <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{110}&space;\widehat{C}(g)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{110}&space;\widehat{C}(g)" title="\widehat{C}(g)" /></a> is the solution of formula (2). During training we iteratively adjust the g vector in the direction of the gradient <a href="https://www.codecogs.com/eqnedit.php?latex=\nabla&space;L" target="_blank"><img src="https://latex.codecogs.com/png.latex?\nabla&space;L" title="\nabla L" /></a>, leading to a g vector, which cellular estimates <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{110}&space;\widehat{C}(g)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{110}&space;\widehat{C}(g)" title="\widehat{C}(g)" /></a> correlate best with the known cellular compositions C. 
