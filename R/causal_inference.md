Causal inference
================
Last updated on 2025-11-24

``` r
library(DiagrammeR) #grViz for causal diagrams
library(dagitty)    #for DAGs
library(ggdag)      #for DAGs
library(tidyverse)


### Set global options
options(digits = 3) # set the default number of digits to 3 


### Rmd settings
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE, fig.path="figures_md/causal_inference/")
```

# Introduction

In this module, we will examine a few case studies to understand the
concept of causal inference. In particular, we focus on understanding
(i) how to correctly estimate the causal effect of the predictor
variable on the response variable, and (ii) whether we can infer the
direction of causality from observational data.

To explore these questions, we will simulate data. The advantage of
simulating data is that we know the true causal relationships between
the variables. This allows us to compare the estimated causal effects
with the true causal effects.

# Case study 1: A confounded relationship

In this case study, we will estimate the correct causal effect of a
predictor variable on a response variable for the following DAG:

``` r
dag <- dagitty('dag {
  X [pos="-1,0"]
  Y [pos="1,0"]
  Z [pos="0, 1"]
  Z -> X
  Z -> Y
}')

ggdag(dag) +
  theme_dag() +
  annotate("text", x = -0.6, y = 0.5, label = "100") +
  annotate("text", x = 0.6, y = 0.5, label = "200")
```

![](figures_md/causal_inference/dag_confound-1.png)<!-- -->

``` r
adjustmentSets(dag, exposure = "X", outcome = "Y")
```

    ## { Z }

In this causal diagram, $Z$ is a confounder that affects both $X$ and
$Y$. The true causal effect of $X$ on $Y$ is 0. To correctly estimate
the causal effect of $X$ on $Y$, we need to adjust for the confounder
$Z$. We will check cases where we do and do not adjust for the
confounder $Z$ and estimate the causal effect of $X$ on $Y$.

Let’s simulate data and estimate the causal effect of $X$ on $Y$.

## Simulate data

``` r
set.seed(123)

### set the effect sizes
beta_zx = 100
beta_zy = 200

### simulate data
n = 10000
z = rbinom(n, 1, 0.6)
x = beta_zx * z + rnorm(n, 500, 100)
y = beta_zy * z + rnorm(n, 500, 100)

data = data.frame(x = x, y = y, z = z)
head(data)
```

    ##     x   y z
    ## 1 551 565 1
    ## 2 613 442 0
    ## 3 485 614 1
    ## 4 648 597 0
    ## 5 592 562 0
    ## 6 634 839 1

## Estimate the confounded causal effect

``` r
m = lm(y ~ x, data = data)
summary(m)
```

    ## 
    ## Call:
    ## lm(formula = y ~ x, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -476.1  -91.4    2.2   92.2  461.9 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 395.7788     6.8014    58.2   <2e-16 ***
    ## x             0.3984     0.0119    33.5   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 133 on 9998 degrees of freedom
    ## Multiple R-squared:  0.101,  Adjusted R-squared:  0.101 
    ## F-statistic: 1.12e+03 on 1 and 9998 DF,  p-value: <2e-16

This confounded model estimates a significant positive association
between $X$ and $Y$, when we set the causal effect of $X$ on $Y$ and
vice verse to be 0. This shows that the confounded model does not
correctly estimate the causal effect of $X$ on $Y$.

## Estimate the causal effect of $Z$ on $X$ and $Y$

``` r
mx = lm(x ~ z, data = data)
summary(mx)
```

    ## 
    ## Call:
    ## lm(formula = x ~ z, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -383.7  -66.8   -0.5   69.7  383.9 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   499.15       1.59   313.3   <2e-16 ***
    ## z             101.73       2.05    49.6   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 100 on 9998 degrees of freedom
    ## Multiple R-squared:  0.198,  Adjusted R-squared:  0.198 
    ## F-statistic: 2.46e+03 on 1 and 9998 DF,  p-value: <2e-16

``` r
my = lm(y ~ z, data = data)
summary(my)
```

    ## 
    ## Call:
    ## lm(formula = y ~ z, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -352.7  -68.6   -0.7   67.6  434.2 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   498.13       1.60   311.4   <2e-16 ***
    ## z             200.11       2.06    97.3   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 101 on 9998 degrees of freedom
    ## Multiple R-squared:  0.486,  Adjusted R-squared:  0.486 
    ## F-statistic: 9.46e+03 on 1 and 9998 DF,  p-value: <2e-16

Both models correctly estimate the causal effect of $Z$ on $X$ and $Y$.

## Estimate the causal effect of $X$ on $Y$ adjusting for $Z$

``` r
m = lm(y ~ x + z, data = data)
summary(m)
```

    ## 
    ## Call:
    ## lm(formula = y ~ x + z, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -353.5  -68.6   -0.5   67.5  434.1 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 492.3224     5.2603   93.59   <2e-16 ***
    ## x             0.0116     0.0100    1.16     0.25    
    ## z           198.9303     2.2966   86.62   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 101 on 9997 degrees of freedom
    ## Multiple R-squared:  0.486,  Adjusted R-squared:  0.486 
    ## F-statistic: 4.73e+03 on 2 and 9997 DF,  p-value: <2e-16

Now, the model with the minimum adjustment set correctly estimates the
causal effect of $X$ on $Y$ to be near 0.

------------------------------------------------------------------------

# Case study 2: Determining the direction of causality

In this case study, we will examine whether we can infer the direction
of causality from observational data. We will consider the following
DAG:

``` r
dag <- dagitty('dag {
  X [pos="-1,0"]
  Y [pos="1,0"]
  X -> Y
}')

ggdag(dag) +
  theme_dag() +
  scale_y_continuous(limits = c(-1, 1)) +
  annotate("text", x = 0, y = 0.1, label = "0.5")
```

![](figures_md/causal_inference/dag_direction-1.png)<!-- -->

In this causal diagram, $X$ causes $Y$. We will estimate the causal
effect of $X$ on $Y$ and $Y$ on $X$ and compare the results.

## Simulate data

``` r
set.seed(123)

### set the effect sizes
beta_xy = 0.5

### simulate data
n = 10000
x = rnorm(n, 500, 100)
y = beta_xy * x + rnorm(n, 500, 100)

data = data.frame(x = x, y = y)
head(data)
```

    ##     x   y
    ## 1 444 959
    ## 2 477 722
    ## 3 656 921
    ## 4 507 697
    ## 5 513 779
    ## 6 672 949

## Estimate the causal effect of $X$ on $Y$

``` r
m = lm(y ~ x, data = data)
summary(m)
```

    ## 
    ## Call:
    ## lm(formula = y ~ x, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -348.3  -66.9   -0.7   68.1  377.0 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  496.071      5.112    97.0   <2e-16 ***
    ## x              0.506      0.010    50.5   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 100 on 9998 degrees of freedom
    ## Multiple R-squared:  0.203,  Adjusted R-squared:  0.203 
    ## F-statistic: 2.55e+03 on 1 and 9998 DF,  p-value: <2e-16

The model correctly estimated the causal effect of $X$ on $Y$ to be 0.5.

Next, we will check the causal effect of $Y$ on $X$, which should be 0.

## Estimate the causal effect of $Y$ on $X$

``` r
m = lm(x ~ y, data = data)
summary(m)
```

    ## 
    ## Call:
    ## lm(formula = x ~ y, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -322.3  -59.9    0.3   60.3  334.9 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 1.99e+02   6.02e+00    33.1   <2e-16 ***
    ## y           4.01e-01   7.95e-03    50.5   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 89.2 on 9998 degrees of freedom
    ## Multiple R-squared:  0.203,  Adjusted R-squared:  0.203 
    ## F-statistic: 2.55e+03 on 1 and 9998 DF,  p-value: <2e-16

The model incorrectly estimated the causal effect of $Y$ on $X$ to be
0.4, when the true causal effect is 0. This shows that we cannot infer
the direction of causality from observational data, at least with
regression models.

We will also check if this holds when we adjust for the confounder $Z$.

------------------------------------------------------------------------

# Case study 3: Determining the direction of causality with a confounder

In this case study, we will examine whether we can infer the direction
of causality from observational data when there is a confounder. We will
consider the following DAG:

``` r
dag <- dagitty('dag {
  X [pos="-1,0"]
  Y [pos="1,0"]
  Z [pos="0, 1"]
  X -> Y
  Z -> X
  Z -> Y
}')

ggdag(dag) +
  theme_dag() +
  annotate("text", x = -0.6, y = 0.5, label = "100") +
  annotate("text", x = 0.6, y = 0.5, label = "200") +
  annotate("text", x = 0, y = 0.05, label = "0.5")
```

![](figures_md/causal_inference/dag_direction_with_confounder-1.png)<!-- -->

## Simulate data

``` r
set.seed(123)

### set the effect sizes
beta_zx = 100
beta_zy = 200
beta_xy = 0.5

### simulate data
n = 10000
z = rbinom(n, 1, 0.6)
x = beta_zx * z + rnorm(n, 500, 100)
y = beta_zy * z + beta_xy * x + rnorm(n, 500, 100)

data = data.frame(x = x, y = y, z = z)
head(data)
```

    ##     x    y z
    ## 1 551  840 1
    ## 2 613  748 0
    ## 3 485  857 1
    ## 4 648  921 0
    ## 5 592  858 0
    ## 6 634 1155 1

## Estimate the causal effect of $X$ on $Y$

``` r
m = lm(y ~ x + z, data = data)
summary(m)
```

    ## 
    ## Call:
    ## lm(formula = y ~ x + z, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -353.5  -68.6   -0.5   67.5  434.1 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  492.322      5.260    93.6   <2e-16 ***
    ## x              0.512      0.010    51.0   <2e-16 ***
    ## z            198.930      2.297    86.6   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 101 on 9997 degrees of freedom
    ## Multiple R-squared:  0.636,  Adjusted R-squared:  0.636 
    ## F-statistic: 8.74e+03 on 2 and 9997 DF,  p-value: <2e-16

The model correctly estimated the causal effect of $X$ on $Y$ to be 0.5
and of $Z$ on $Y$ to be 200.

Next, we will check the causal effect of $Y$ on $X$, which should be 0.

## Estimate the causal effect of $Y$ on $X$

``` r
m = lm(x ~ z, data = data)
summary(m)
```

    ## 
    ## Call:
    ## lm(formula = x ~ z, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -383.7  -66.8   -0.5   69.7  383.9 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   499.15       1.59   313.3   <2e-16 ***
    ## z             101.73       2.05    49.6   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 100 on 9998 degrees of freedom
    ## Multiple R-squared:  0.198,  Adjusted R-squared:  0.198 
    ## F-statistic: 2.46e+03 on 1 and 9998 DF,  p-value: <2e-16

``` r
m = lm(x ~ y, data = data)
summary(m)
```

    ## 
    ## Call:
    ## lm(formula = x ~ y, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -315.8  -60.8    0.6   60.4  333.5 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 1.97e+02   4.90e+00    40.2   <2e-16 ***
    ## y           4.04e-01   5.35e-03    75.5   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 89.3 on 9998 degrees of freedom
    ## Multiple R-squared:  0.363,  Adjusted R-squared:  0.363 
    ## F-statistic: 5.7e+03 on 1 and 9998 DF,  p-value: <2e-16

``` r
m = lm(x ~ y + z, data = data)
summary(m)
```

    ## 
    ## Call:
    ## lm(formula = x ~ y + z, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -315.9  -60.8    0.6   60.4  333.5 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 1.98e+02   6.08e+00   32.52   <2e-16 ***
    ## y           4.03e-01   7.91e-03   50.96   <2e-16 ***
    ## z           5.75e-01   2.70e+00    0.21     0.83    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 89.3 on 9997 degrees of freedom
    ## Multiple R-squared:  0.363,  Adjusted R-squared:  0.363 
    ## F-statistic: 2.85e+03 on 2 and 9997 DF,  p-value: <2e-16

Although the first model without y correctly estimated the causal effect
of $Z$ on $X$, the second model with y incorrectly estimated the causal
effect of $Y$ on $X$ to be 0.4, when the true causal effect is 0, as
well as of $Z$ on $X$, which should be 100. This shows that we cannot
infer the direction of causality from observational data, even when
adjusting for the confounder $Z$. In addition, it shows that
conditioning on a collider can introduce bias in the causal effect
estimates.

<br>

------------------------------------------------------------------------

# Case study 4: Effect of using the same variable for calculating rate

In this case study, we explore whether calculating the rates using the
same variable introduces spurious correlation.

## Simulate data

``` r
set.seed(123)

### simulate data
n = 10000
x = rnorm(n, 500, 100)
y = rnorm(n, 500, 100)
z = rnorm(n, 500, 100)
x_rate = x / z
y_rate = y / z

data = data.frame(x = x, y = y, z = z, 
                  x_rate = x_rate, y_rate = y_rate)
head(data)
```

    ##     x   y   z x_rate y_rate
    ## 1 444 737 416  1.066   1.77
    ## 2 477 483 478  0.998   1.01
    ## 3 656 593 290  2.264   2.05
    ## 4 507 443 333  1.522   1.33
    ## 5 513 523 390  1.315   1.34
    ## 6 672 613 333  2.014   1.84

## Check the correlation

``` r
cor(data$x, data$y)
```

    ## [1] 0.00602

``` r
cor(data$x_rate, data$y_rate)
```

    ## [1] 0.557

``` r
cor(data$x, data$y_rate)
```

    ## [1] -0.0135

This shows that calculating rates based on the same variable leads to
spurious correlation.
