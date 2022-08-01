--- 
title: "Methods in (Skene & Kenward, 2010)"
author: "Dylan Dijk"
date: "2022-08-01"
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








```r
servr::daemon_stop(1)
```


# Introduction

This tutorial aims to show how methods described in the two papers by @skene_analysis_2010_I^,^^[@skene_analysis_2010_II] can be applied with mixed models in [R](#R) or in [SAS](#SAS).

In particular we will focus on: 

  * Modified Box Correction
  * Scheffe contrast
  * *Can add more as I go*

***

## SAS

The analysis for these two papers were carried out using SAS. 
  
## R

There are multiple packages in R that allow you to fit mixed models. The most popular are the **nlme** and **lme4** package. Some of the differences between the two packages are described on page 4 of the **lme4** package documentation (https://cran.r-project.org/web/packages/lme4/lme4.pdf). 

The main difference is that the **nlme** package allows more flexibility in choosing complex variance-covariance structures. With the **lme4** package only allowing us to select the correlation structure between individuals via the random effects (G side).


