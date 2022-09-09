
# Scheffé’s method {#Scheffe-method}





## Code
For the Scheffé test of the individual contrasts we use values computed in the modified Box correction test. Therefore we use the `Box_correction_modified_nlme()` function from previously within another function. 






```r
Scheffe = function(model, model_r, data, Sigma = NULL, a){
  
  x = Box_correction_modified_nlme(model, model_r, data, Sigma, print = FALSE)
  

C_a = sum((x$fixed_effects_estimate + (c(0,rep(x$fixed_effects_estimate[1],6)))) * a)
sigma_2 = (t(x$Y) %*% x$A %*% x$Y)/(x$n-x$r)

n_i = 3

S_Ca = sqrt(sigma_2 * sum(a^2/n_i))
S_alpha_a = S_Ca * sqrt((x$c * x$Lambda) * qf(0.95,x$c,x$v2_mod))

print(C_a)
print(paste0("(",C_a - S_alpha_a,", ", C_a + S_alpha_a, ")"))
  
}
```




## Example

We have already shown the calculation of the modified box corrected statistic for the guinea pig example in Section \@ref(Guinea-pigs). We now extend this to calculate the confidence intervals given in Table X of paper II.



```r
Guinea_pigs_data = read.csv("https://raw.githubusercontent.com/DylanDijk/Simon-Skene-Mixed-Models-Tutorial/main/Data_Images_Figures/brammer.csv")
```

We also need to make sure that the concentration variable is stored as a factor variable:

```r
Guinea_pigs_data$conc.f = factor(Guinea_pigs_data$conc.f)
```




```r
library(nlme)

model_AP1 = gls(AP1~conc.f,data=Guinea_pigs_data)
model_AP1_r = gls(AP1~1,data=Guinea_pigs_data)

model_AP2 = gls(AP2~conc.f,data=Guinea_pigs_data)
model_AP2_r = gls(AP2~1,data=Guinea_pigs_data)
```


```r
sigma_compound_1 = cov(rbind(Guinea_pigs_data[1:7,1],Guinea_pigs_data[8:14,1],Guinea_pigs_data[15:21,1]))
sigma_compound_2 = cov(rbind(Guinea_pigs_data[1:7,2],Guinea_pigs_data[8:14,2],Guinea_pigs_data[15:21,2]))
```

Next we create the blocked version of these covariance matrices, a block for each of the three subjects.


```r
m = 3 #no. subjects

sigma_compound_1_block = diag(m) %x% sigma_compound_1
sigma_compound_2_block = diag(m) %x% sigma_compound_2

```



```r
Scheffe(model = model_AP1, model_r = model_AP1_r, data = Guinea_pigs_data, Sigma = sigma_compound_1_block, a = c(-1,1,0,0,0,0,0))
```
For compound 1, we get these values:

```
[1] 1.333333
[1] "(-3.07944941452269, 5.74611608118969)"
```



