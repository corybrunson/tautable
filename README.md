tautable
=========

Scripts to produce tables of critical values for Kendall's tau-a

## Install

```r
require(devtools)
if(!require(tautable)) {
    install_github('corybrunson/tautable')
    if(!require(tautable)) {
        stop('package is not available')
    }
}
```
