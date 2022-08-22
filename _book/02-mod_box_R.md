# (PART\*) R {.unnumbered}


# Modified Box Correction

## nlme 



Here we provide a function that calculates the modified box correction from mixed models constructed using the **nlme** package.

The function takes as input two **nlme** models, the input variable `model` is the full model then `model_r` should be the model with the variables removed that we want to test for significance in the model.



```r

Box_correction_modified_nlme = function(model, model_r, data, Sigma = NULL){
  
  
  
  # Design matrix for full model
  X <- model.matrix(model,data = data)
  
  ind = which(is.na(match(rownames(data),rownames(X))))

  # Response variable vector
  Y = getResponse(model)

  
  # Design matrix for model with removed variables
  X_r <- model.matrix(model_r,data = data)
  # Number of observations
  n = nrow(X)
  # Number of variables
  parm = ncol(X)
  # Number of variables being removed
  c = ncol(X)-ncol(X_r)


  # Have added option for Sigma to be chosen manually, this is needed for the Guinea pig example.
  # By default the function takes Sigma to be NULL, so will be NULL unless manually defined.
  
  if(is.null(Sigma)){
  
      # Need to find good way of extracting covariance matrix from model
        # Put in if statement to deal with models with no correlated errors
      if(is.null(model$modelStruct$corStruct)){
        
        Sigma = diag((model$sigma)^2, nrow = n, ncol = n)
      
      } else {
        # repeats = as.numeric(levels(model$groups))
        
        # Creating block matrix for model with correlated errors
        Sigma = getVarCov(model, type = "marginal")
        Sigma = diag(nlevels(model$groups)) %x% Sigma
          if(!(length(ind) == 0)){
            Sigma = Sigma[-ind,-ind]
          }
        }
     
    }
     
  
  A = diag(n) - ( X %*% solve( t(X) %*% X ) %*% t(X) )
  B = ( X %*% solve( t(X) %*% X ) %*% t(X) )  -  ( X_r %*% solve( t(X_r) %*% X_r ) %*% t(X_r) ) 
  
  trace = function(x){
    sum(diag(x))
  }
  
  V = ( (trace((B %*% Sigma) %*% (B %*% Sigma)) / trace(B %*% Sigma)^2)
        + (trace((A %*% Sigma) %*% (A %*% Sigma)) / trace(A %*% Sigma)^2) )   
  
  
  v_2 = (c*(4*V + 1) - 2)/((c*V) - 1)
  
  lambda = ((n - parm)/c) * ((v_2 - 2)/v_2) * ((trace(B %*% Sigma) / trace(A %*% Sigma)))
  
  F_ =  ((n - parm)/c) * ((t(Y) %*% B %*% Y)/(t(Y) %*% A %*% Y))
  
  F_mod = F_/lambda
  
  p_value = pf(df1 = c, df2 = v_2, q =  F_mod, lower.tail = F)
  
  values = c(lambda, F_, F_mod, c, v_2, p_value)
  names(values) = c("Lambda", "F_","F_MOD", "v1_MOD", "v2_mod", "prob_F_mod")
  values = round(values, digits = 3)
  
  return(as.table(values))
  
}

```


In your R session run this code so that you then have the function in your R environment.

Next we give some examples of using the function to recreate results presented in the second paper.

## Examples

In this section we look at recreating the examples given in section 6 of the second paper. 

### Cardiac enzyme

This example is covered in Section 6.1 of the second  paper^[@skene_analysis_2010_II]. Section 6.1 of the paper applies the modified Box correction twice to this dataset, once with the original versiona and after including artificial dropouts to the dataset.

Importing the dataset from GitHub repository, and converting a selection of variables to factors.


```r
Cardiac_enzyme_data = read.csv("https://raw.githubusercontent.com/DylanDijk/Simon-Skene-Mixed-Models-Tutorial/main/Data_Images_Figures/The%20Cardiac%20Enzyme%20Data%20-%20Reduced%20Data%20set.csv")
Cardiac_enzyme_data$dog = as.factor(Cardiac_enzyme_data$dog)
Cardiac_enzyme_data$trt = as.factor(Cardiac_enzyme_data$trt)
Cardiac_enzyme_data$time = as.factor(Cardiac_enzyme_data$time)
```

Here we are creating the **nlme** models. In this examples I have set `glsControl(tolerance = 10^-8)`, this is too match the convergence criteria used by **nlme** to the default convergence criteria used by PROC MIXED in SAS.


```r
library(nlme)
glsControl(tolerance = 10^-8)

model = gls(model = atp ~  trt + time + trt:time,
            data = Cardiac_enzyme_data,
            correlation = corSymm(form =~ as.numeric(time)|dog),
            weights = varIdent(form =~ 1|time))


model_r = gls(model = atp ~ trt + time,
              data = Cardiac_enzyme_data,
              correlation = corSymm(form =~ as.numeric(time)|dog),
              weights = varIdent(form =~ 1|time))
```



We can compare the results from the `Box_correction_modified_nlme()` function results given in the upper panel of Table VII in the second paper. 


```r
Box_correction_modified_nlme(data = Cardiac_enzyme_data, model = model, model_r = model_r)
```


 
In addition to computing the p-value using the modified Box corrected F statistic, we can use the **lmerTest** package to compute some of the other p-values in Table VII.  

The other p-values are computed using Wald tests with the Kenward-Roger adjustment.

Here we 

```r
library(lmerTest)
m1 = lmerTest::lmer(atp ~  trt + time + trt:time + (1|dog),
              data = Cardiac_enzyme_data)

m1 = lme4::lmer(atp ~  trt + time + trt:time + (1|dog),
              data = Cardiac_enzyme_data)

anova(m1, ddf = "Kenward-Roger")
```
The anova method for the S3 method for class 'lmerModLmerTest' uses the KRmodcomp function from the pbkrtest package. And as an alternative to the previous method we can use the KRmodcomp function directly.


```r
library(pbkrtest)
m0 = lme4::lmer(atp ~ trt + time + (1|dog),
                data = Cardiac_enzyme_data)
KRmodcomp(m1,m0)

```




### Cardiac enzyme - Artifical missing data

The results in the lower panel of Table VII are computed by applying the same method to the same dataset but now with artificial missing data introduced.

The `Box_correction_modified_nlme()` function works with missing data. 

Below we give the R code to compute these results.

We first introduce the missing data points to the `Cardiac_enzyme_data` used in the previous section.


```r
Cardiac_enzyme_data_miss = Cardiac_enzyme_data
Cardiac_enzyme_data_miss$atp[which(Cardiac_enzyme_data_miss$dog == 4)][7:9] = NA
```


```r
library(nlme)
glsControl(tolerance = 10^-8)


model = gls(model = atp ~  trt + time + trt:time,
            data = Cardiac_enzyme_data_miss,
            correlation = corSymm(form =~ as.numeric(time)|dog),
            weights = varIdent(form =~ 1|time), na.action = na.omit)

model_r = gls(model = atp ~ trt + time,
              data = Cardiac_enzyme_data_miss,
              correlation = corSymm(form =~ as.numeric(time)|dog),
              weights = varIdent(form =~ 1|time), na.action = na.omit)


Box_correction_modified_nlme(model = model, model_r = model_r, data = Cardiac_enzyme_data_miss)
```

### Guinea pigs {#Guinea-pigs}

Importing the dataset.


```r
Guinea_pigs_data = read.csv("Data_Images_Figures/brammer.csv")

Guinea_pigs_data$conc.f = factor(Guinea_pigs_data$conc.f)

```

For this example @skene_analysis_2010_II used the sample covariance matrix as an estimator of the covariance matrix ($\Sigma$). Therefore the sample covariance estimates are computed and then are used as inputs for `Sigma` in the `Box_correction_modified_nlme()` function.


Calculating the sample covariance matrices.

```r
sigma_compound_1 = cov(rbind(Guinea_pigs_data[1:7,1],Guinea_pigs_data[8:14,1],Guinea_pigs_data[15:21,1]))
sigma_compound_2 = cov(rbind(Guinea_pigs_data[1:7,2],Guinea_pigs_data[8:14,2],Guinea_pigs_data[15:21,2]))
```



```r
m = 3 #no. subjects

sigma_compound_1_block = diag(m) %x% sigma_compound_1
sigma_compound_2_block = diag(m) %x% sigma_compound_2

```


Creating the models.

```r
library(nlme)

model_AP1 = gls(AP1~conc.f,data=Guinea_pigs_data)
model_AP1_r = gls(AP1~1,data=Guinea_pigs_data)

model_AP2 = gls(AP2~conc.f,data=Guinea_pigs_data)
model_AP2_r = gls(AP2~1,data=Guinea_pigs_data)
```



Running the `Box_correction_modified_nlme()` function. 

  * For compound 1.

```r
Box_correction_modified_nlme(model = model_AP1, model_r = model_AP1_r, data = Guinea_pigs_data, Sigma = sigma_compound_1_block)
```
  * For compound 2.

```r
Box_correction_modified_nlme(model = model_AP2, model_r = model_AP2_r, data = Guinea_pigs_data, Sigma = sigma_compound_2_block)
```
  


