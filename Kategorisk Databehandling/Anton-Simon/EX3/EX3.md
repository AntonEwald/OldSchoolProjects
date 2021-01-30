Computer assignment 2: Logistic regression
================
Simon Melamed & Anton Holm
11/21/2019

### Summary

In this first assignment we are given data regarding people with
periodontitis collected from a group of adult patients from a dentist
clinic. The patients were compared to a control group with normal gums
who also visited the clinic. At the same time, patients were asked about
whether they floss or not. The data are gathered in the table below.

|                       | Periodontitis: Yes | Periodontitis: No | Sum |
| :-------------------: | :----------------: | :---------------: | :-: |
| **Dental floss: Yes** |         22         |        75         | 97  |
| **Dental floss: No**  |        148         |        265        | 413 |
|        **Sum**        |        170         |        340        | 510 |

# Exercise 2.1

### Task 1

We now fit the logistic regression model of the probability for
periodontitis in the population as a function of dental floss. This
entails that we have \(Y\) as the response variable and \(X\) the
covariate. We assert the logistic regression model:

\[
\text{logit}(p_x) = \beta_0 + \beta_1 x 
\]

where,

\[
p_x = P(Y=1 | X=x) = \frac{exp(\beta_0 + \beta_1 x)}{1+exp(\beta_0 + \beta_1 x)} \\
x = 0: \text{Dental floss is not used}\\
x = 1: \text{Dental floss is used} \\
\]

were \(Y\) and \(X\) being categorial variables and the support of \(Y\)
being \(0,1\) where \(Y=0\) means the patient do not have periodontitis
and \(Y=1\) the patient does. Below follows an R output for the logistic
regression model.

    ## 
    ## Call:
    ## glm(formula = y ~ x, family = binomial(link = logit), data = data21a, 
    ##     weights = n)
    ## 
    ## Deviance Residuals: 
    ##       1        2        3        4  
    ## -15.335   17.429   -6.212    8.080  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -0.5825     0.1026  -5.677 1.37e-08 ***
    ## x            -0.6439     0.2633  -2.446   0.0145 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 649.24  on 3  degrees of freedom
    ## Residual deviance: 642.80  on 2  degrees of freedom
    ## AIC: 646.8
    ## 
    ## Number of Fisher Scoring iterations: 5

First of all we can see that both the intercept \(\beta_0\) and the
effect parameter \(\beta_1\) are significant for the model on the
significance level \(0.01\).

This type of study is a retrospective case-control study. Since we
decided how many patients to have in both groups (having or not having
periodontitis) and \(X\) only having two levels we get binomial
distributed columns. This leads to us not being able to estimate
\(\beta_0\) which we will talk more about further down. \(\beta_1\) on
the other hand is not a problem. Assume that we, in a more general
scenario of having \(X\) as an interval scale predictor variable,
increase \(x\) to \(x+1\). We then get the logarithm odds ratio
\(\theta\) as:

\[
log(\theta) = \log\frac{(p_{x+1})(1-p_{x+1})}{(p_x)(1-p_x)} = logit(p_{x+1})-[~logit(p_x)~] = \beta_1
\]

So the effect parameter \(\beta_1\) is the logarithm of the odds ratio
\(\theta\) between \(X=x+1\) and \(X=x\). In our case with \(x\) only
taking the values \(0\) and \(1\) we have that \(\beta_1\) is the odds
ratio of the odds for \(X=1\) and the odds for \(X=0\). Having an
estimation of \(\beta_1\) as \(\hat{\beta}_1 = -0.6439\), which implies
that the estimated odds of getting periodontitis when using dental floss
is \(e^{\hat{\beta}_1}=0.52524\) the odds of getting periodontitis when
not using dental floss.

Let us now consider \(\beta_0=logit(p_0)\). To get an estimation of
\(\beta_0\) we can, when no cell observations are small, use
\(\hat{\beta}_0=y_0/n_0\). Here, \(y_i\) is the observed value of a
binomial distributed variable \(Y_i\) with success probobility \(p_i\)
and \(n_i\) is the ammount of tries at the level \(X=x_i\). For the case
when some cells don’t have many observation we merge the observations in
that cell with another cell. The problem here lies in the fact that the
data is sampled such that we have binomial distributed columns with a
fixed number of patients with periodontitis and a fixed number of
case-control patients. Hence the amount of people using dental floss or
not (\(n_1\) & \(n_0\)) are random and change depending on what patients
we chose for both groups.

Another way to put it is that \(\beta_0\) depends on the sampling
probabilities \(q_0\) and \(q_1\). From Bayes Formula we get, after some
calculation, that

\[
P(Y=1~|~~sampled,x)=...=\frac{exp(\beta_0^\star + \beta_1 x)}{1+exp(\beta_0^\star+\beta_1x)}=p_x(\beta_0^\star,\beta_1)
\]

where \(\beta_0^\star=\beta_0+log(p_1/p_0)\) and \(q_i\) is:
\(q_i = P(\text{sampled}~|~Y=i, X=x)\)

which requires us to know the total number of controls and cases in the
whole population. If we had this information we could estimate \(q_i\)
as \(q_i=n_i/N_i\) where \(n_i\) is the number of controls (\(i=0\)) and
number of cases (\(i=1\)) in our study, \(N_0\) the total number of
controls in the entire population and \(N_1\) the total number of cases
in the entire population. Luckily, this has no impact on \(\beta_1\) so
we can still draw conclusions from our data.

### Task 2

In the second task, we regard \(X\) as the response variable instead and
\(Y\) as the predictor. This gives us a logistic regression model of the
probability for using dental floss depending on if you have
periodontisis or not, which entails:

\[
logit (p_y) = \gamma_0 + \gamma_1 y \\
y=0~~ \text{ if you don't have periodontisis}\\
y=1~~ \text{ if you have periodontisis}\\
p_y = P(X=1~|~ Y=y)
\]

The summary of information regarding this logistic regression model can
be seen below:

    ## 
    ## Call:
    ## glm(formula = x ~ y, family = binomial(link = logit), weights = n)
    ## 
    ## Deviance Residuals: 
    ##       1        2        3        4  
    ## -11.493   -6.405   15.057    9.485  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -1.2622     0.1308  -9.651   <2e-16 ***
    ## y            -0.6439     0.2633  -2.446   0.0145 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 496.24  on 3  degrees of freedom
    ## Residual deviance: 489.79  on 2  degrees of freedom
    ## AIC: 493.79
    ## 
    ## Number of Fisher Scoring iterations: 5

Just as in the previous excercise, both variables are significant on a
\(0.01\) significance level meaning that \(\gamma_1 \neq 0\) on a
\(99 \%\) confidence level and \(\gamma_0 \neq 0\) actually on an even
higher confidence level.

In this model we do not have any problems estimating the intercept. We
get \(\gamma_0\) as:

\[
\gamma_0 = logit(p_{y=0}) = log\frac{p_{y=0}}{1-p_{y=0}}
\]

with,

\[
p_{y=0}=P(X=1~|~Y=0)\approx\frac{n_{12}}{n_{+2}}=\frac{75}{340}
\]

which in the end gives us an estimation of \(\gamma_0\) as:

\[
\gamma_0 \approx \log( \frac{75}{340}\cdot (1-\frac{75}{340})^{-1} )= \log(0.2830189)= -1.262242.
\]

\(\gamma_1\) is this time, in the same way as in the previous excercise,
the odds ratio between \(Y=y+1\) and \(Y=y\). The estimated odds-ratio
for a \(2\)x\(2\) table doesn’t regard the conditioning of \(X\) and
\(Y\) but only takes into consideration the cell observations which is
why the estimation of the effect parameter in both logistic regression
models are the same.

# Exercise 2.2

### Summary

For this exercise we were given data from a prospective study where 160
mice were exposed to different doses of bensoapyren, which is a
polycyclic aromatic hydrocarbon which is suspected to enhance the risk
of devoloping lung tumors. The results from the study is shown in the
table below, the x column is the number of lung tumors detected in that
trial, whereas the logdos column is the log of the dose of besoapyren
and the n column is the total amount of mice in each trial.

| logdos | x  | n  |
| :----: | :-: | :-: |
| \-7.6  | 1  | 18 |
| \-6.22 | 2  | 19 |
| \-4.6  | 4  | 28 |
|  \-3   | 9  | 32 |
| \-1.39 | 12 | 28 |
|  0.92  | 32 | 40 |

### Task 1

For the second execise, we calculated the estimates of the risk to
develop a lung tumor for each dose of bensoapyren, \(p_i =P(Y=1|X=x_i)\)
where \(i \in \{1,2,3,4,5,6\}\) for each dose of bensoapyren. We can
estimate these probabilities by dividing the amount that devoloped tumor
by the total amount that were given the dose,
\(p_i = \frac{n_{i1}}{n_{i+}}\). Using this method, the retrieve the
estimates to be,
\[p_1 = \frac{1}{18}, p_2 = \frac{2}{19}, p_3 = \frac{4}{28}, p_4 = \frac{9}{32}, p_5 = \frac{12}{28}, p_6 = \frac{32}{40}.\]
A plot of these estimates against the logarithm of the dose can be found
in Plot 1.

If we would like to calculate the estimated odds, \(O_i\) for each dose
we simply calculate \(\frac{p_{i}}{1-p_i}\) för each dose
\(x_i, i \in \{1,2,3,4,5,6\}\). The result of this method yield,
\[O_1 = \frac{\frac{1}{18}}{1-\frac{1}{18}},O_2 = \frac{\frac{2}{19}}{1-\frac{2}{19}},O_3 = \frac{\frac{4}{28}}{1-\frac{4}{28}},O_4 = \frac{\frac{9}{32}}{1-\frac{9}{32}},O_5 = \frac{\frac{12}{28}}{1-\frac{12}{28}},O_6 = \frac{\frac{32}{40}}{1-\frac{32}{40}}.\]
We’ve plotted the logarithm of these odds against the log of the dose of
bensoapyren, see Plot 2. In this case, when we apply a logistic
regression, we have that
\[\text{logit}(p_{x_{i}})= \beta_0 +\beta_1x_i\] So, for a logistic
model to be appropriate for this dataset, the slope of the Plot 2 should
resemble a linear line, which it does. Therefore, we can apply a
logistic regression model to the data.

![](EX3_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->![](EX3_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

### Task 2

Now, when we have concluded through plotting that logistic regression
seems to be suitable for the dataset we can start to construct a model
and calculate estimates. We assert the logistic regression model,

\[
\text{logit}(p_x) = \beta_0 + \beta_1 x 
\] where,

\[
p_x = P(Y=1 | X=x) = \frac{e^{\beta_0 + \beta_1 x}}{1+e^{\beta_0 + \beta_1 x}} \\
x = x_i: \text{Log dose of bensoapyren}.\\
\] For this model, the risk of getting a lung tumor is the response
variable and log of the dose of bensoapyren is the covariate. If we plug
in this model to R and analyze the output we can see that we get the
coefficients \[\beta_0 =  0.68702 , \beta_1 = 0.52037,\] with p-values
below the significance level of \(\alpha = 0.01.\) With this, we can
state, with \(99 \%\) statistical confidence, that exposure of
bensoapyren have an impact on the risk of devoloping lung tumors.

The interpretation of the \(\beta_1\) coefficient is, that for every
unit increase of the dose in the log scale, the log of the odds of
getting a tumor increases with \(0.52037\) units. However, in this
scale, an increase of \(0.52037\) units is not very intuitive and
difficult to interpret. Another way of thinking of it, is that the log
of the odds ratio of getting a tumor equal
\(log(\theta)=\beta_1 = 0.52037\), which implies
\(\theta = e^{0.52037}= 1.68265\). This interprets as, when the log of
the dose increases with 1 unit, the odds of getting a lung tumor
increases with approximately \(68 \%\). So, the steeper the slope, or
the bigger \(\beta_1\), the greater impact the dose of bensoapyren have
on the odds of getting a lung tumor.

We can plug in our estimates into our risk function,

\[
p_x = P(Y=1 | X=x) = \frac{e^{  0.68702+0.52037 x}}{1+e^{  0.68702+0.52037 x}} \\
x = x_i: \text{Log dose of bensoapyren}.\\
\] With this, we can calculate the risk of getting a lung tumor for
different doses of bensoapyren. For example, the median of the log dose,
which for this study is \(x = -3.8\) which implies a dose of
\(e^{-3.8}\approx 0.02237077\) , which asserts,

\[
p_x = P(Y=1 | X=-3.8) = \frac{e^{  0.68702+0.52037 \cdot (-3.8)}}{1+e^{  0.68702+0.52037 \cdot (-3.8)}} =0.2157875.
\]

So the risk of getting a lung tumor when the log dose of besnoapyren is
\(-3.8\), is approximately \(22 \%\).

The interpretation of \(\beta_0=0.68702\) is that, when the log dose of
bensoapyren equal to zero, \(log(Dose)=0\) which is equivalent to the
dose being equal to one, \(Dose=e^0=1\), we have that the risk of
getting a tumor is approximately \(67 \%\),

\[
p_x = P(Y=1 | X=-3.8) = \frac{e^{ 0.68702}}{1+e^{0.68702}} = 0.6653037.
\]

### Task 3

Below is the covariance matrix for the logistic regression model
previously discussed, calculated using R’s built in function.

    ##             (Intercept)      logdos
    ## (Intercept)  0.06633368 0.014358467
    ## logdos       0.01435847 0.007229089

From this we can calculate the estimatedmcorrelation between \(\beta_0\)
and \(\beta_1\), which is given by known formula,
\[Cor(\beta_0,\beta_1)=\frac{Cov(\beta_0,\beta_1)}{\sqrt{Var(\beta_0) \cdot Var(\beta_1)}}=0.655691.\]
This implies a positive correlation between the intercept \(\beta_0\)
and the effect parameter \(\beta_1\)

Now that we have interpreted and calculated the estimates, we can
continue to find \(95\%\) confidence interval for them. This is easily
done by using the information gathered from the covariance matrix. We
find in the output that the standard deviation (square root of the
variance) for \(\beta_0\) is \(0.25755\) and so \(95\%\) Wald-confidence
for \(\beta_0\) interval is given by

\[(\hat{\beta}_0-1.96\cdot SE_{\beta_0},\hat{\beta}_0+1.96\cdot SE_{\beta_0})=(0.182222, 1.191818),\]

where \(1.96\) is the \(97.5\%\) and \(2.5\%\)-quantile for the standard
normal distribution. Similarily, it is found in the output that the
standard deviation for \(\beta_1\) is \(0.08502\) and so a a \(95%\)
Wald-confidence for \(\beta_1\) interval is given by
\[(\hat{\beta}_1-1.96\cdot SE_{\beta_1},\hat{\beta}_1+1.96\cdot SE_{\beta_1})=(0.3537308, 0.6870092)\].
From this we can calculate the conficence interval for the effect of one
unit increase in log dose of bensoapyren on the odds of getting a tumor,
\((e^{0.3537308},e^{0.6870092})=(1.424372,1.987762).\) With this we can
state, on a \(95\%\) confidence level, that an increase of one unit of
log dose, increase the odds of getting a lung tumor with at least
\(42 \%\) and with at most \(99\%\).

To find a \(95 \%\) confidence interval for the tumor risk at dose
\(0.25\), meaning \(\log(Dose)= -1.39)\), we need to find a confidence
interval for \(\beta_0+\beta_1x\). But then we need to calculate the
standard deviation for \(\beta_0+\beta_1x\). This is done by taking the
square root of the variance, \(Var(\beta_0+\beta_1x)\). The variance is
retrieved by the following computation,

\[Var(\beta_0+\beta_1x)=Var(\beta_0)+x^2Var(\beta_1)+2xCov(\beta_0,\beta_1).\]
The variance for \(\beta_0\) and \(\beta_0\), along with the covariance,
\(Cov(\beta_0,\beta_1)\) is found in the covariance matrix above. With
this we calculate the variance, when \(x=-1.39\) to be,
\[Var(\beta_0+\beta_1x)=0.06633368+(-1.39)^2 \cdot 0.007229089 + 2(-1.39) \cdot  0.01435847=0.04038446,\]
and so the standard deviation becomes \(\sqrt{0.04038446}=0.2009588\).
The \(95\%\) Wald confidence interval for \(\beta_0+\beta_1\cdot x\) is
then given by,
\[(\beta_0+\beta_1\cdot x + 1.96 \cdot SE_{\beta_0+\beta_1\cdot x},\beta_0+\beta_1\cdot x - 1.96 \cdot SE_{\beta_0+\beta_1\cdot x}) =  (-0.4301735,0.3575849).\]
To find the confidence interval for the risk of developing a lung tumor
when log dose of bensoapyren is \(-1.39\), we simply compute,
\[(\frac{e^{-0.4301735}}{1+e^{-0.4301735}},(\frac{e^{0.3575849}}{1+e^{0.3575849}})\approx (0.394, 0.588).\]
With this, we can state with \(95 \%\) confidence level that the risk of
developing a lung tumor when log dose of bensoapyren is \(-1.39\), lies
approximately between \(39 \%\) and \(59 \%\).

### Task 4

Since the Wald confidence inteval for \(\beta_1\) does not contain zero
we can suspect that, when contructing a hypothesis test for
\(H_0: \beta_1 = 0\) and \(H_1: \beta_1 \neq 0\), we should be able to
conclude that we should believe \(H_1\) in favor of \(H_0\). To properly
research this, we compute a Wald-test and a likelihood ratio test.

The test statistic for the Wald test can be construced as,
\[z_w^2= \frac{\hat{\beta_1}^2}{Var(\hat{\beta_1})} \stackrel{H_0}{\sim} \chi^2(1).\]
Note that we’ve used the square of the Wald-statistic since these values
are easily collected using the information gathered from the summary
output, \[z_w^2= \frac{0.52037^2}{0.007229089}=37.45768.\] However, one
could have calculated the square root of this result as well. The
Wald-statistic chosen is for large enough samples converging towards a
chi-squared distribution so we compare the statistic with the \(95 \%\)
quantile of a chi-squared distribution. Since the \(95 \%\) quantile for
the chi-squared distribution with one degree of freedom is \(3.84\), and
our test statistic is much larger than the quantile, we can reject the
null hypothesis in favor of the alternative on the \(95\%\) confidence
level using the wald test. This statistic can also be found in the
summary output as the square of the z-value for \(\beta_1\).

The test statistic for the likelihood ratio test is constructed as,

\[-2(L_p(\beta_0(\beta_1=0),0)- L(\hat{\beta_0}, \hat{\beta_1})) \stackrel{H_0}{\sim} \chi^2(1),\]
where \(L_p(\beta_0(\beta_1=0),0)\) is the profile log-likelihood under
the null hypothesis that \(\beta_1 =0\). The profile log-likelihood is
calculated by maximizing the log-likelihood while \(\beta_1\) is fixed,
for this case is equal to zero, so \(\beta_0(\beta_1=0)\) is the value
of \(\beta_0\) that maximize the log-likelihood function when
\(\beta_1=0\). And \(L(\hat{\beta_0}, \hat{\beta_1})\) is simply the
maximum log-likelihood with our estimates as parameters. The likelihood
ratio statistic also converges towards a chi-squared distribution, in
this case with one degree of freedom. Using the code given in the
instructions, we can see below that the likelihood ratio test statistic
is \(55.291\), which also is larger than the \(95 \%\) quantile for the
chi-squared distribution with one degree of freedom (\(55.291 > 3.84\)).
So the interpretation of this test is also that, we can, on a \(95 \%\)
confidence level, reject the null hypothesis that \(H_0: \beta_1 = 0\)
in favor of \(H_1: \beta_1 \neq 0\).

    ## Single term deletions
    ## 
    ## Model:
    ## x/n ~ logdos
    ##        Df Deviance    AIC    LRT Pr(>Chi)    
    ## <none>       1.242 24.024                    
    ## logdos  1   56.532 77.315 55.291 1.04e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Task 5

In this task we aim to expand our understanding of how the sample size
affect the variance of the estimated parameters. Since we have no
expression explicitly for the variances, we cannot confirm it with the
use of a formula. However, in this task we’ll research how the sample
sizes affect the variances when we multiply the table used in Exercise 2
with 10, 100 and 1000 respectively.

    ##                                   Beta 0       Beta 1
    ## Original sample size        6.633368e-02 7.229089e-03
    ## Original sample size x 10   6.633369e-03 7.229093e-04
    ## Original sample size x 100  6.633369e-04 7.229093e-05
    ## Original sample size x 1000 6.633369e-05 7.229093e-06

It becomes clear that for each time the counts in the table is
multiplied by 10, the variance for each parameter is scaled down with 10
(devided by 10). This seems reasonable, as a larger sample size should
intuitively result in a better estimation and a smaller confidence
interval for each parameters.
