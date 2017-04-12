# PGCA: An Algorithm to Link Protein Groups Created from MS/MS Data

## Installation
The latest version of the R package can be installed directly from this
repository. First the dependencies and tools have to be installed into R.
```r
# Install dependencies
install.packages(c("knitr", "devtools"))
```

When the dependencies are fulfilled, the `devtools` package can be used
to install the `pgca`:
```r
# Install pgca
library(devtools)
install_github("gcohenfr/pgca", build_vignettes = TRUE)
```

## Usage
To get a quick introduction how to use the package, simply look at the
vignette:
```r
vignette("intro", package = "pgca")
```