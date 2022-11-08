# lemad
### Install 

To install the R package LEMAD, make sure you have devtools installed and then type from R (Rstudio):

```
library(devtools)
remotes::install_github("leonelhalsina/lemad")
```
It is likely that the latest Rtools needs to be installed.
We recommend to have the latest R (it works perfectly with R 4.2.1).
Also, package rgal might cause issues.


### Using Lemad

We have prepared a vignette (a sort of manual with chunks of code) that can
be called once you load the library:

```
library(lemad)
browseVignettes("lemad")
```
### News

1- On November 8th we have released a version where c++ code (from Thijs Janzen, University of Groningen) is used to reduce computation time. It is a shame that we could not include this nice feature at the time of publishing the paper three months ago. In the paper we mentioned that a 66-species tree and 3 areas (which means it is total of 7 areas because all combinations) takes 10 minutes to complete analysis. With this c++ improvement, that time is shorted. It now takes 5 seconds.


### Cite

When using lemad, please cite:

Herrera-Alsina L, Algar A. C., Lancaster L.T., Ornelas J. F., Bocedi G, Papadopulos A. S. T., Gubry-Rangin C., Osborne O. G., Mynard P., Sudiana I. M., P Lupiyaningdyah P., Juliandi B. & Travis J.M.J. (2022). The Missing Link in Biogeographic Reconstruction: Accounting for Lineage Extinction Rewrites History. Journal of Biogeography https://doi.org/10.1111/jbi.14489
