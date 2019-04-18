# DTD
Digital Tissue Deconvolution (DTD) reconstructs the cellular composition of a tissue from its bulk expression profile. 
In order to increase deconvolution accuracy, DTD adapts the deconvolution model to the tissue scenario via loss-function learning. 
Training is performed on 'in-silicio' training mixtures, for which the cellular composition are known.  
As input, DTD requires a labelled expression matrix. 
The package includes functions to generate training and test mixtures, train the model, and assess its deconvolution capability via visualizations.  
In Goertler et al. 2018 "Loss-function Learning for Digital Tissue Deconvolution" the theory has been published 

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
