--- 
title: "Methods in (Skene & Kenward, 2010)"
author: "Dylan Dijk"
date: "2022-09-09"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
# url: your book url like https://bookdown.org/yihui/bookdown
# cover-image: path to the social sharing image like images/cover.jpg
description: |
  This is a minimal example of using the bookdown package to write a book.
  The HTML output format for this example is bookdown::bs4_book,
  set in the _output.yml file.
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---

<!-- I have followed the bookdown documentation found at https://bookdown.org/yihui/bookdown/ to make this.  -->

<!-- This section from this tutorial has been useful https://happygitwithr.com/bookdown-cheat-sheet.html -->



<!-- When I want to link to GitHub repository look at section 3.1.2.7 in Bookdown documentation -->




```r
servr::daemon_stop(1)
```


# Introduction

This tutorial aims to show how methods described in the [two papers](#References) by @skene_analysis_2010_I^,^^[@skene_analysis_2010_II] (paper I and paper II) can be applied in [R](#R-Intro).

<!-- For the rest of this tutorial we will refer to these two papers as paper I and paper II respectively. -->

In these papers, it is assumed that the data can be represented by a multivariate Gaussian linear model. The model has the following form:

\begin{equation}
  y_i \sim N(X_i \beta;\Sigma_i ), \quad i =1, \dots,n
  (\#eq:general-model)
\end{equation}


<!-- where $y_i$ $(T_i \times 1)$ is the response vector from the $i$th subject of $n$ subjects. $X_i$ $(T_i \times p)$ the design matrix and $\Sigma_i$ $(T_i \times T_i)$ the covariance matrix for subject $i$.  -->

where $y_i$ $(T_i \times 1)$ is the response vector, $X_i$ $(T_i \times p)$ the design matrix and $\Sigma_i$ $(T_i \times T_i)$ the covariance matrix for the $i$th subject of $n$ subjects.



Papers I and II cover various methods that can be used for analysis of repeated measurements under this framework, and below is a quick overview of some of these methods.

<!-- Equivalently this can be written as: -->

<!-- \[ y \sim N(X \beta;\Sigma ) \] -->

<!-- where $y = (y_1,\dots, y_n)^T$, $X = (X_1,\dots, X_n)^T$ and $\Sigma = \text{block-diagonal}\{\Sigma_i\}$ -->

## Methods we will focus on: 

### Adjusted Sandwich Estimator {#Overview-Adjusted-Sandwich-Estimator}

<!-- In this set up the estimate of the covariance matrix influences both the estimate of $\beta$ and then the estimation of the variance of $\beta$. -->

The sandwich estimator, shown in paper I, avoids the estimate of $\Sigma$ in the calculation of the estimate of $\beta$. Then uses an empirical estimate of $\Sigma$ to give an estimate of the variance of $\hat\beta$. 

However, the tests based on statistics using this estimator are unreliable. In Section 3 of paper I an $F$-test using a Wald statistic which takes account of the variability in the sandwich estimator is shown.



### Modified Box Correction  {#Overview-Modified-Box-Correction}
The Box correction is a correction of the one way ANOVA $F$-statistic which follows an $F$ distribution under the assumptions of independent equally variable Gaussian errors. Therefore we depart from the general model given by \@ref(eq:general-model) and we assume $\Sigma = \text{diag}(\sigma^2)$.

The Box correction $\psi^{-1}$$F$ then follows an $F$ distribution under the more general model, where the Gaussian errors no longer are restricted to being independent and equally variable.

In the **modified Box correction** the ratio of quadratic forms are treated as a scaled $F$-statistic instead of a ratio of independent chi-squared distributions. The modified Box correction is preferred to Box???s original statistic which is conservative and hence less powerful. 

In order to compute the modified Box correction we require a consistent estimator of the covariance matrix $\Sigma$. There are many options for this, for example, you could use the **OLS** covariance estimate. However, in most software that fit mixed models the **REML** estimator of the covariance matrix is computed and is therefore often more practical to use this as an estimate. The ideal case is to compute the unstructured REML estimate, but when this is not possible we can compute an estimate of $\Sigma$ with the most complex structure the data will support.

In addition, it is also possible to use an empirical estimator such as the sample  covariance matrix, an example of this is shown in Section 6.2 of paper II. We recreate the results for this example with R in subection \@ref(Guinea-pigs).

    
    
    
### Scheff?? contrast {#Overview-Scheffe-contrast}
  
Alongside the modified Box correction to test for the significance of covariates on the response variable, paper II also covers how Scheff?????s method can be adapted to use the modified Box corrected statistic.

Scheff?????s method can be used to test for any possible contrasts between means, and in paper II Scheff?????s method is used to give confidence intervals for individual contrasts in the guinea pig papillary muscle example (see Section 6.2 of paper II). We recreate the results of this example with R in Section \@ref(Scheffe-method).  

***

## SAS {#SAS-Intro}

The analysis for these two papers were carried out using SAS. In particular, the procedure **PROC MIXED** is used to calculate an estimator of the covariance matrix for the modified Box correction. SAS provides [documentation](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_mixed_toc.htm "SAS PROC MIXED documentation") for **PROC MIXED**, and in the penultimate paragraph within the [overview section](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_mixed_overview.htm) it describes how this SAS procedure fits the structure you select to the data by using the method of restricted maximum likelihood (REML).

In this tutorial we focus on providing functions in R to compute the methods described in the papers, and we only present the SAS code used by Simon Skene to calculate the results for example 6.1 in paper II. 

***  
  
## R {#R-Intro}

For all three methods listed above ([Sandwich Estimator](#Sandwich-Estimator), [Modified Box Correction](#Modified-Box-Correction) and [Scheff?? contrasts](#Scheffe-method)) we provide R functions that compute them. The modified Box correction requires an estimate for $\Sigma$ and in paper II a REML estimate is used. In order to compute this estimate we utilise R packages that fit mixed models.

There are multiple packages that fit mixed models in R, with the most popular being the **nlme** and **lme4** packages. Both of these packages use REML by default to fit the mixed model, therefore we can use them to extract an estimate of the covariance matrix to compute the modified Box correction.

Some of the differences between the two packages are described on page 4 of the **lme4** package [documentation](https://cran.r-project.org/web/packages/lme4/lme4.pdf "lme4 R package documentation"). The main difference is that the **nlme** package allows more flexibility in choosing complex variance-covariance structures. On the other hand, the **lme4** package only allows us to induce a correlation structure between individuals by including random effects in the model.


