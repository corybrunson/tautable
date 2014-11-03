tautable
=========

Scripts to produce tables of critical values for Kendall's tau-a

## Install

```r
require(devtools)
if(!require(tautable)) {
    install_github('corybrunson/tautable')
    if(!require(taut able)) {
        stop('package is not available')
    }
}
```
