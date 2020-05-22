Optimization of threshold parameters for GFA replicates

### Installation
To install in R:

```
install.packages("devtools")
devtools::install_github("kforthman/optmThrGFA")
library(optmThrGFA)
```

To install on the server:

```
# first, log onto the server through the terminal
ssh submit0
# (type your password and hit 'enter')
module load R
R
install.packages('devtools')
# You may need to select a mirror. I tend to use the one for Kansas. Type the number and press 'enter'.
# Once you select your mirror, it may take a couple of minutes to look like it's doing anything,
# so just be patient.
devtools::install_github("kforthman/optmThrGFA")
# This step is only to check the package installed correctly:
library(optmThrGFA)
```

If you have trouble installing any dependencies on the server, you may need to install an earlier version of the package that works with R 3.5.1. You may do this using the R code below
```
library(devtools)
devtools::install_version("package_name", version = "package_version")
```
