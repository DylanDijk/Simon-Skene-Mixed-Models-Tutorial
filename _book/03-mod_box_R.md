


# Modified Box Correction {#Modified-Box-Correction}

## nlme 



The code below gives an R function that returns the modified Box correction and the corresponding p-value. The function uses mixed model objects from the **nlme** package to define the mean model, and the structure of the covariance matrix that will be used. The preferred choice is to use an unstructured estimate, but when this is not possible we can compute an estimate of $\Sigma$ with the most complex structure the data will support. 


The function requires two **nlme** models as input, the `model_r` variable should be the formula with the variables removed that you want to test significance for in the full model defined by the `model` variable. 

The default estimate of $\Sigma$ is the REML estimate extracted from the model object defined by the `model` input variable. However, if you do want to use an alternative covariance estimate you can select it using the `Sigma` variable. For example in Section \@ref(Guinea-pigs) we use the sample covariance matrix.



```r
library(nlme)

Box_correction_modified_nlme = function(model, model_r, data, Sigma = NULL){
  
  rownames(data) = c()
  # Design matrix for full model
  X <- model.matrix(model,data = data)
  
  # Gives the rownumbers of the dataset that have missing values
  ind = which(is.na(match(rownames(data),rownames(X))))
  
  Y = getResponse(model)
  
  # Design matrix for model with removed variables
  X_r <- model.matrix(model_r,data = data)
  # Number of observations
  n = nrow(X)
  # Number of variables
  parm = ncol(X)
  # Number of variables being removed
  c = ncol(X)-ncol(X_r)
  
  if(is.null(Sigma)){
    
    # Put in if statement to deal with models with no correlated errors
    if(is.null(model$modelStruct$corStruct)){
      
      Sigma = diag((model$sigma)^2, nrow = n, ncol = n)
      
    } else {
      # repeats = as.numeric(levels(model$groups))
      
      # Creating block matrix for model with correlated errors
      # type = "marginal" extracts the covariance matrix when 
      # not conditioning on the random effects. Therefore extracts 
      # the full covariance matrix G side and R side.
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
  
  values = c(lambda, F_, c, v_2, F_mod, p_value)
  names(values) = c("Lambda", "F_", "v1_MOD", "v2_mod","F_MOD", "prob_F_mod")
  values = round(values, digits = 3)
  
  
  fixed_effects_estimate = solve(t(X)%*%X)%*%t(X)%*%Y
  
  return(list("Lambda" = lambda, "F_" = F_[1], "c" = c, "v2_mod" = v_2, "F_MOD" = F_mod[1],       
              "prob_F_mod" = p_value[1], "A" = A, "B" = B, "n" = n, "r" = parm, "Y" = Y,
              "fixed_effects_estimate" = fixed_effects_estimate))
  
}
```


In your R session copy and run this code so that it is available in your R environment. Next we give some examples of applying the function to recreate results presented in the paper II.

## Examples

### Cardiac enzyme

In this example we look at the same dataset and mean model that was used in the example covered in Section 6.1 of paper II. We apply the `Box_correction_modified_nlme()` function to test for the inclusion of the interaction term between treatment and time in the model.

We apply the function with the original dataset and then to the dataset with artificial dropouts, to recreate the results of Table VII.

First of all we need to import the dataset. The code below imports the dataset from a GitHub repository, and converts a selection of the variables in the dataset to factors.


```r
Cardiac_enzyme_data = read.csv("https://raw.githubusercontent.com/DylanDijk/Simon-Skene-Mixed-Models-Tutorial/main/Data_Images_Figures/The%20Cardiac%20Enzyme%20Data%20-%20Reduced%20Data%20set.csv")
Cardiac_enzyme_data$dog = as.factor(Cardiac_enzyme_data$dog)
Cardiac_enzyme_data$trt = as.factor(Cardiac_enzyme_data$trt)
Cardiac_enzyme_data$time = as.factor(Cardiac_enzyme_data$time)
```

The code below creates the **nlme** models, a good resource documenting how to fit different mixed models with **nlme** is [Fitting linear mixed models in R](http://staff.pubhealth.ku.dk/~jufo/courses/rm2018/nlmePackage.pdf) by Brice Ozenne.

I have also set `glsControl(tolerance = 10^-8)`, this is too match the convergence criteria used by **nlme** to the default convergence criteria used by PROC MIXED in SAS.


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

 <!-- model_r_2 = gls(model = atp ~ trt + time, -->
 <!--               data = Cardiac_enzyme_data, -->
 <!--               correlation = corSymm(form =~ as.numeric(time)|dog)) -->

We can compare the results from the `Box_correction_modified_nlme()` function with the results given in the upper panel of Table VII in paper II. 


```r
cardiac_mod_box = 
Box_correction_modified_nlme(data = Cardiac_enzyme_data, model = model, model_r = model_r)
```


```r
cardiac_mod_box$F_MOD
[1] 2.520353
cardiac_mod_box$prob_F_mod
[1] 0.07743055
```

<!-- # The part below shows how to calculate other p-values in table VII -->
<!-- In addition to computing the p-value using the modified Box corrected F statistic, we can use the **lmerTest** package to compute some of the other p-values in Table VII.   -->

<!-- The other p-values are computed using Wald tests with the Kenward-Roger adjustment. -->

<!-- Here we  -->
<!-- ```{r} -->
<!-- library(lmerTest) -->
<!-- m1 = lmerTest::lmer(atp ~  trt + time + trt:time + (1|dog), -->
<!--               data = Cardiac_enzyme_data) -->

<!-- m1 = lme4::lmer(atp ~  trt + time + trt:time + (1|dog), -->
<!--               data = Cardiac_enzyme_data) -->

<!-- anova(m1, ddf = "Kenward-Roger") -->
<!-- ``` -->
<!-- The anova method for the S3 method for class 'lmerModLmerTest' uses the KRmodcomp function from the pbkrtest package. And as an alternative to the previous method we can use the KRmodcomp function directly. -->

<!-- ```{r} -->
<!-- library(pbkrtest) -->
<!-- m0 = lme4::lmer(atp ~ trt + time + (1|dog), -->
<!--                 data = Cardiac_enzyme_data) -->
<!-- KRmodcomp(m1,m0) -->

<!-- ``` -->




### Cardiac enzyme - Artifical missing data

The results in the lower panel of Table VII are computed by applying the same method to the same dataset but now with artificial missing data introduced. The `Box_correction_modified_nlme()` function works with missing data, and we run it with the modified dataset in the code below.

We first introduce the missing data points to the `Cardiac_enzyme_data` dataset:


```r
Cardiac_enzyme_data_miss = Cardiac_enzyme_data
Cardiac_enzyme_data_miss$atp[which(Cardiac_enzyme_data_miss$dog == 4)][7:9] = NA
```

The code below creates the models with the modified dataset, it is identical to the model construction in the previous example but with the new dataset.


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
```


We get the same resuls as the final row of Table VII in paper II:

```r
cardiac_miss_mod_box = 
Box_correction_modified_nlme(model = model, model_r = model_r, data = Cardiac_enzyme_data_miss)
```


```r
cardiac_miss_mod_box$F_MOD
[1] 2.84006
cardiac_miss_mod_box$prob_F_mod
[1] 0.05912901
```



### Guinea pigs {#Guinea-pigs}

In this example we use a dataset that shows the amplitude of action potential that was measured from pieces of heart dissected from guinea pigs. The measurements were taken after being exposed to two different compounds, therefore we have two response variables AP1 and AP2.

Each piece of tissue had a control measure where no compound was given and then six increasing concentrations of the compound. This was carried out on three guinea pigs' hearts, so we have seven repeated measurements from just three participants. 

The code below imports the dataset and converts the concentration variable (`conc.f`) and subject variable (`gpig.f`) to factors:


```r
Guinea_pigs_data = read.csv("https://raw.githubusercontent.com/DylanDijk/Simon-Skene-Mixed-Models-Tutorial/main/Data_Images_Figures/brammer.csv")

Guinea_pigs_data$conc.f = factor(Guinea_pigs_data$conc.f)
Guinea_pigs_data$gpig.f = factor(Guinea_pigs_data$gpig.f)
```


In this example the experiments are treated as simple repeated measures designs with concentration as the time variable and tissue as the subject. 

Due to the small number of samples, the estimation method for the unstructured covariance models did not converge. Therefore, for this example Skene and Kenward used the sample covariance matrix as an estimator of the covariance matrix ($\Sigma$).

The code below computes the sample covariance estimates for the seven repeated measurements, which are then used as inputs for `Sigma` in the `Box_correction_modified_nlme()` function.


```r
sigma_compound_1 = cov(rbind(Guinea_pigs_data[1:7,"AP1"],Guinea_pigs_data[8:14,"AP1"],Guinea_pigs_data[15:21,"AP1"]))
sigma_compound_2 = cov(rbind(Guinea_pigs_data[1:7,"AP2"],Guinea_pigs_data[8:14,"AP2"],Guinea_pigs_data[15:21,"AP2"]))
```

Next we create the blocked version of these covariance matrices, a block for each of the subjects. To create this covariance matrix we need to know the number of unique subjects in the dataset.


```r
number_of_subjects = nlevels(Guinea_pigs_data$gpig.f) #number of subjects
 
sigma_compound_1_block = diag(number_of_subjects) %x% sigma_compound_1
sigma_compound_2_block = diag(number_of_subjects) %x% sigma_compound_2
```


Next we construct the models, we make a model for each of the response variables with concentration (`conc.f`) as the only covariate.

```r
library(nlme)

model_AP1 = gls(AP1~conc.f,data=Guinea_pigs_data)
model_AP1_r = gls(AP1~1,data=Guinea_pigs_data)

model_AP2 = gls(AP2~conc.f,data=Guinea_pigs_data)
model_AP2_r = gls(AP2~1,data=Guinea_pigs_data)
```


Below we run the `Box_correction_modified_nlme()` function, using the model objects we have created above as inputs. 

  * For compound 1.

```r
Box_mod_1 =
Box_correction_modified_nlme(model = model_AP1, model_r = model_AP1_r, data = Guinea_pigs_data, Sigma = sigma_compound_1_block)
```
This gives the following output:

```r
Box_mod_1$F_MOD
[1] 4.497613
Box_mod_1$prob_F_mod
[1] 0.05698665
```
  * For compound 2.

```r
Box_mod_2 = Box_correction_modified_nlme(model = model_AP2, model_r = model_AP2_r, data = Guinea_pigs_data, Sigma = sigma_compound_2_block)
```
For compound 2, we get these values:

```r
Box_mod_2$F_MOD
[1] 20.84697
Box_mod_2$prob_F_mod
[1] 0.001995718
```



