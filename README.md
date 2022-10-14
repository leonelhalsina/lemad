# lemad
**Install** 

To install the R package LEMAD, make sure you have devtools installed and then type from R (Rstudio):

```
library(devtools)
remotes::install_github("leonelhalsina/lemad")
```
It is likely that the latest Rtools needs to be installed.
We recommend to have the latest R (it works perfectly with R 4.2.1).
Also, package rgal might cause issues.

**Using Lemad**

We have prepared a vignette (a sort of manual with chunks of code) that can
be called once you load the library:

```
library(lemad)
browseVignettes("lemad")
```
**Citing Lemad**

When using lemad, please cite:

Herrera-Alsina L, Algar A. C., Lancaster L.T., Ornelas J. F., Bocedi G, Papadopulos A. S. T., Gubry-Rangin C., Osborne O. G., Mynard P., Sudiana I. M., P Lupiyaningdyah P., Juliandi B. & Travis J.M.J. (2022). The Missing Link in Biogeographic Reconstruction: Accounting for Lineage Extinction Rewrites History. Journal of Biogeography https://doi.org/10.1111/jbi.14489
