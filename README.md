This repository contains an R package for performing "Variance Adaptive Shrinkage"

To install the vashr package first you need to install devtools
```
install.packages("devtools")
library(devtools)
install_github("mengyin/vashr")
```

## Running Variance Adaptive Shrinkage


The main function in the ashr package is `vash`. To get minimal help:
```
> library(vashr)
> ?vash
```

## More background

The vashr ("Variances Adaptive SHrinkage") package aims to provide a flexible approach to derive ``shrinkage-based" estimates for unknown variances $s=(s_1,\dots,s_G)$, given only the estimated standard errors ($\hat{s}=(\hat{s}_1,\dots,\hat{s}_G)$). Our key idea is to combine information across the multiple measurements $g=1,\dots,G$ to improve the estimation accuracy for each individual $s_g$.  

Here we propose an adaptive approach to shrink the variances, such that the appropriate amount of shrinkage is determined from the data. Our method is based on the assumption that the distribution of the variances (or, alternatively, the precisions) is unimodal, which allows more flexibility than the widely used package ["limma"](https://bioconductor.org/packages/release/bioc/html/limma.html). We use a mixture of (possibly a large number of) inverse-gamma distributions to flexibly model this unimodal distribution, and provide simple computational procedures to fit this Empirical Bayes model by maximum likelihood of the mixture proportions.

Our procedure provides a posterior distribution on each variance or precision, as well as point estimates (posterior mean). The methods are an analogue of the ``adaptive shrinkage" methods for mean parameters introduced in [Stephens (2016),"False Discovery Rates: A New Deal"](http://biorxiv.org/content/early/2016/01/29/038216). 
