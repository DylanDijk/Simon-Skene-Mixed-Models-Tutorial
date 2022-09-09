
# Scheffé’s method {#Scheffe-method}




## Code

The code below gives an R function that uses Scheffé's method with the modified Box corrected statistic to calculate the confidence intervals for individual contrasts within a categorical effect.
To use this `Scheffe()` function we require the `Box_correction_modified_nlme()` function to be in our R environment. 

The `Scheffe()` function requires the input variables `model` and `model_r`, in the same way it is used for the `Box_correction_modified_nlme()` function. However, as we are looking at contrasts within an individual categorical variable `model_r` should only have a single categorical variable removed relative to the formula selected by `model`. For instance, in the [example below](#Scheffe-Guinea-pig-example) `model_r` has the categorical variable `conc.f` removed.

This function has an input variable `a` which is used to select the contrasts that we want to calculate confidence intervals for. Equation (22) in paper II shows how this $a$ vector is used to define the contrasts.




```r
Scheffe = function(model, model_r, data, Sigma = NULL, a, subjects){
  
  x = Box_correction_modified_nlme(model, model_r, data, Sigma)
  
C_a = sum((x$fixed_effects_estimate + (c(0,rep(x$fixed_effects_estimate[1],6)))) * a)
sigma_2 = (t(x$Y) %*% x$A %*% x$Y)/(x$n-x$r)

cat_var_removed = setdiff(attr(model$terms, "term.labels"), attr(model_r$terms, "term.labels"))

n_i = vector()
for(i in 1:nlevels(data[[cat_var_removed]])){
  n_i[i] = sum(!is.na(data[data[[cat_var_removed]] == levels(data[[cat_var_removed]])[i],cat_var_removed]))
}

S_Ca = sqrt(sigma_2 * sum(a^2/n_i))
S_alpha_a = S_Ca * sqrt((x$c * x$Lambda) * qf(0.95,x$c,x$v2_mod))

return(list("C_a" = C_a, "Conf_interval" = paste0("(",C_a - S_alpha_a,", ", C_a + S_alpha_a, ")")))
  
}
```


## Example {#Scheffe-Guinea-pig-example}

We have already shown the calculation of the modified box corrected statistic for the guinea pig example in subsection \@ref(Guinea-pigs). We now extend this to calculate the confidence intervals given in Table X of paper II.

The setup before we run the `Scheffe()` function to get the confidence intervals is identical to the setup in subsection \@ref(Guinea-pigs).




```r
Guinea_pigs_data = read.csv("https://raw.githubusercontent.com/DylanDijk/Simon-Skene-Mixed-Models-Tutorial/main/Data_Images_Figures/brammer.csv")

Guinea_pigs_data$conc.f = factor(Guinea_pigs_data$conc.f)
Guinea_pigs_data$gpig.f = factor(Guinea_pigs_data$gpig.f)
```




```r
library(nlme)

model_AP1 = gls(AP1~conc.f,data=Guinea_pigs_data)
model_AP1_r = gls(AP1~1,data=Guinea_pigs_data)

model_AP2 = gls(AP2~conc.f,data=Guinea_pigs_data)
model_AP2_r = gls(AP2~1,data=Guinea_pigs_data)
```


```r
sigma_compound_1 = cov(rbind(Guinea_pigs_data[1:7,"AP1"],Guinea_pigs_data[8:14,"AP1"],Guinea_pigs_data[15:21,"AP1"]))
sigma_compound_2 = cov(rbind(Guinea_pigs_data[1:7,"AP2"],Guinea_pigs_data[8:14,"AP2"],Guinea_pigs_data[15:21,"AP2"]))
```


```r
number_of_subjects = nlevels(Guinea_pigs_data$gpig.f) #number of subjects

sigma_compound_1_block = diag(number_of_subjects) %x% sigma_compound_1
sigma_compound_2_block = diag(number_of_subjects) %x% sigma_compound_2
```

The choice of `a` and the contrasts we want to look at is based on the ordering of the levels of the categorical variable. 
Therefore if we want to compute the confidence interval for the mean difference between the control measure and the first concentration, we set `a = c(-1,1,0,0,0,0,0)`.

For compound 1, we get these values:


```r
Scheffe_comp_1 = 
Scheffe(model = model_AP1, model_r = model_AP1_r, data = Guinea_pigs_data, Sigma = sigma_compound_1_block, a = c(-1,1,0,0,0,0,0), subjects = "gpig.f")
```


```r
Scheffe_comp_1$Conf_interval
[1] "(-3.07944941452269, 5.74611608118969)"
```


<!-- ```{r, echo=FALSE, results='markup'} -->
<!-- Scheffe_comp_2 =  -->
<!-- Scheffe(model = model_AP2, model_r = model_AP2_r, data = Guinea_pigs_data, Sigma = sigma_compound_2_block, a = c(-1,1,0,0,0,0,0)) -->
<!-- ```    -->



