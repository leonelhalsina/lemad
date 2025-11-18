# lemad
### Version
This is LEMAD v2.1.1

### Install 

To install the R package LEMAD, make sure you have devtools installed. There is an issue with the latest RcppParallel package, so 
I recomend to install an older, more stable version of it (5.1.9). During the installation of lemad, you might be
asked to update RcppRparallel from 5.1.9 to the latest, please skip that update. In R (Rstudio) type:

```
library(devtools)
remotes::install_version("RcppParallel", version = "5.1.9")
remotes::install_github("leonelhalsina/lemad")
```
### If installation gives you a hard time:
It is likely that the latest Rtools needs to be installed.
We recommend to have the latest R (it works perfectly with R 4.2.2 and 4.5.2, tested in both
Windows 10 and Unix systems).
Also, package rgal might cause issues.

In some Unix systems, it might need the following line when installing the vignette:
pandoc::pandoc_install()

For MAC users, it is being reported that installing gfortran 11.2 could
be useful.

DO NOT hesitate in contact me if none of this works: leonelhalsina@gmail.com


### Using Lemad

We have prepared a vignette (a sort of manual with chunks of code) that can
be called once you load the library:

```
library(lemad)
devtools::install(build_vignettes = TRUE)
# or perhaps you need to build the vignette since package installation, so do:
remotes::install_github("leonelhalsina/lemad",build_vignettes = TRUE)
browseVignettes("lemad")
```
### Wiki
Have a look at our Wiki to see news and a record of important upgrades to LEMAD. For instance, on November 2022 LEMAD became super fast: https://github.com/leonelhalsina/lemad/wiki

### Cite

When using lemad, please cite:

Herrera-Alsina L, Algar A. C., Lancaster L.T., Ornelas J. F., Bocedi G, Papadopulos A. S. T., Gubry-Rangin C., Osborne O. G., Mynard P., Sudiana I. M., P Lupiyaningdyah P., Juliandi B. & Travis J.M.J. (2022). The Missing Link in Biogeographic Reconstruction: Accounting for Lineage Extinction Rewrites History. Journal of Biogeography https://doi.org/10.1111/jbi.14489
