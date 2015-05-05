pkgname <- "JSM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('JSM')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("JSM-package")
### * JSM-package

flush(stderr()); flush(stdout())

### Name: JSM-package
### Title: What the package does (short line) Joint modeling
### Aliases: JSM-package JSM
### Keywords: Joint Modeling

### ** Examples

1+3



cleanEx()
nameEx("jmodelMult")
### * jmodelMult

flush(stderr()); flush(stdout())

### Name: jmodelMult
### Title: Joint Modeling Main Function with NMRE (nonparametric
###   Multiplicative random effects)
### Aliases: jmodelMult

### ** Examples

1 + 3



cleanEx()
nameEx("jmodelTM")
### * jmodelTM

flush(stderr()); flush(stdout())

### Name: jmodelTM
### Title: Joint Modeling Main Function with LME (linear mixed effects)
### Aliases: jmodelTM

### ** Examples

1 + 3



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
