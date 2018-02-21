# dcm

Implements inverse probability weighting (IPW) with the logit discrete choice nonresponse model (Tchetgen Tchetgen et. al, 2017).

# How to Install

To install this package, use devtools:

```r
devtools::install_github("bluosun/dcm")
```
(requires R version >= 3.4.1) 

and then load the package:

```r
library("dcm")
```

# Overview
This function implements IPW using LDCM weights as described in Tchetgen Tchetgen et. al (2017). The implementation is based on 
a default linear main effects model of the observed variables in the *r*<sup>th</sup> missing data pattern for the log ratio 
of probabilities for observing the *r*<sup>th</sup> missing data pattern versus the complete data.

# Usage
 ```r
 dcm(data_frame, regr_formula, regr_fam)
 ```
# Arguments
* data_frame: A data frame containing the variables in the regression model.
* regr_formula: An object of class "formula" describing the regression model to be fitted
* regr_fam: A description of the error distribution and link function to be used in the regression model. This can be a character string (e.g. "binomial" or "gaussian") naming a family function, same as the input to glm().

# Return values
A "dcm" object containing the following items:
* CC: An object of class "glm" from complete-case regression analysis.
* IPW: The corresponding coefficient point estimate, bootstrap standard error, confidence interval and p-value from IPW analysis.
* DAT: A data frame containing the variables in the regression model and the missing data pattern indicator R.

# Example
The 'airquality' dataset is included in R base and contains nonmonotone missing values in the variables Ozone (ppb) and Solar.R (lang). 
```r
out <- dcm(airquality, Temp~Ozone+Solar.R+Wind, gaussian)
```

# Reference

Tchetgen Tchetgen, E., Wang, L. and Sun, B. (2017). <a href="http://www3.stat.sinica.edu.tw/ss_newpaper/SS-2016-0325_na.pdf"> Discrete Choice Models for Nonmonotone Nonignorable Missing Data: Identification and Inference.</a> Statistica Sinica (doi: 10.5705/ss.202016.0325).

