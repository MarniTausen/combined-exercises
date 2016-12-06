Project exercises 3
================
Marni Tausen
2016-12-06

Class project: Bayesian linear regresssion
------------------------------------------

Implement a constructor for a `blm` class. One approach, taken from the textbook, is implementing an `update` function and a `blm` function:

``` r
update <- function(model, prior, ...) { ... }
blm <- function(model, ...) {
    # some code here...
    prior <- make_a_prior_distribution_somehow()
    posterior <- update(model, prior, ...)
    # some code that returns an object here...
}
```

To get this version of `blm` to work you need to get the prior in a form you can pass along to `update` but if you did the exercises earlier you should already have a function that does this (although you might want to create a class for these distributions and return them as such so you can manipulate them through an interface if you want to take it a bit further).

``` r
make_prior <- function(model, ...){

    parameters <- list(...)
    if(!is.null(parameters$alpha)) { alpha <- parameters$alpha }
    else { alpha <- 1 }

    model <- model.matrix(model)
    list(mu=vector("numeric", ncol(model)), Sigma=diag(1/alpha, nrow=ncol(model)))
}

update <- function(model, prior, alpha, beta, ...) {
    data <- model.frame(model)
    mx <- model.matrix(model, ...)
    Sxy <- solve(diag(alpha, nrow=ncol(mx)) + beta * t(mx) %*% mx)
    mxy <- beta * Sxy %*% t(mx) %*% data$y
    list(mu=mxy, Sigma=Sxy)
}


blm <- function(model, alpha, beta, ...) {

    prior <- make_prior(model, ...)
    posterior <- update(model, prior, alpha, beta, ...)

    obj <- list(coefficients=t(posterior$mu),
                variances=posterior$Sigma,
                model=model.frame(model),
                prior=prior,
                posterior=posterior,
                terms=model,
                alpha=alpha,
                beta=beta)
    class(obj) <- "blm"
    obj
}
```

``` r
x <- 1:100
y <- rnorm(100, x)
z <- x/y

fit <- blm(y ~ x + z, 1, 1)
```

### Model methods

There are some polymorphic functions that are generally provided by classes that represent fitted models. Not all models implement all of them, but the more you implement, the more existing code can manipulate your new class; another reason for providing interfaces to objects through functions only.

Below is a list of functions that I think your `blm` class should implement. The functions are listed in alphabetical order, but many of them are easier to implement by using one or more of the others. So read through the list before you start programming. If you think that one function can be implemented simpler by calling one of the others, then implement it that way.

In all cases, read the R documentation for the generic function first. You need the documentation to implement the right interface for each function anyway so you might at least read the whole thing. The description in this note is just an overview of what the functions should do.

#### coefficients

This function should return fitted parameters of the model. It is not entirely straightforward to interpret what that means with our Bayesian models where a fitted model is a distribution and not a single point parameter. We could let the function return the fitted distribution, but the way this function is typically used that would make it useless for existing code to access the fitted parameters for this model as a drop in replacement for the corresponding parameters from a `lm` model, for example. Instead, it is probably better to return the point estimates of the parameters which would be the mean of the posterior you compute when fitting.

Return the result as a numeric vector with the parameters named. That would fit what you get from `lm`.

``` r
coefficients.blm <- function(x) x$coefficients

coefficients(fit)
```

    ##      (Intercept)        x         z
    ## [1,]    1.561226 0.998953 -1.577903

#### confint

The function `confint` gives you confidence intervals for the fitted parameters. Here we have the same issue as with `coefficients`: we infer an entire distribution and not a parameter (and in any case, our parameters do not have confidence intervals; they have a joint distribution). Nevertheless, we can compute the analogue to confidence intervals from the distribution we have inferred.

If our posterior is distributed as **w** ∼ *N*(**m**, **S**) then component *i* of the weight vector is distributed as *w*<sub>*i*</sub> ∼ *N*(*m*<sub>*i*</sub>, **S**<sub>*i*, *i*</sub>). From this, and the desired fraction of density you want, you can pull out the thresholds that match the quantiles you need.

You take the `level` parameter of the function and get the threshold quantiles by exploiting that a normal distribution is symmetric. So you want the quantiles to be `c(level/2, 1-level/2)`. From that, you can get the thresholds using the function `qnorm`.

``` r
#    parm: a specification of which parameters are to be given
#          confidence intervals, either a vector of numbers or a vector
#          of names.  If missing, all parameters are considered.
confint.blm <- function(object, parm, level = 0.95, ...){

    fit <- fitted(object)
}
```

#### deviance

This function just computes the sum of squared distances from the predicted response variables to the observed. This should be easy enough to compute if you could get the squared distances, or even if you only had the distances and had to square them yourself. Perhaps there is a function that gives you that?

``` r
deviance.blm <- function(object, ...){
    (object$model[,1]-fitted(object, ...)$mean)^2
}

deviance(fit)
```

    ##   [1] 2.278178e-01 8.602888e-02 3.882410e-02 9.819216e-03 4.759771e-02
    ##   [6] 4.378543e+00 1.364795e+00 4.538412e-02 1.006626e-01 2.536898e+00
    ##  [11] 4.017819e-01 5.882770e-01 6.048392e-01 5.938171e-02 1.295352e+00
    ##  [16] 7.338530e-03 7.311772e-01 7.279517e-03 2.768422e-03 2.194194e+00
    ##  [21] 3.078178e-01 1.677081e-01 3.404886e-02 1.943493e+00 9.775070e-03
    ##  [26] 1.043138e-01 8.449061e-01 2.842467e-01 1.566926e+00 2.516129e-01
    ##  [31] 2.859780e+00 1.472052e+00 3.333890e-02 1.862363e+00 1.471825e+00
    ##  [36] 4.802963e-03 5.611307e-02 7.524098e-01 3.521478e-02 3.202501e-04
    ##  [41] 1.435698e+00 1.696999e-01 6.453958e-02 1.170179e+00 1.491195e+00
    ##  [46] 5.190291e-01 2.291314e-02 2.577161e+00 4.540375e-02 3.361012e+00
    ##  [51] 7.616651e-02 9.588254e-01 4.752506e-01 1.823319e+00 9.855274e-01
    ##  [56] 3.432933e-04 9.134839e-02 4.405019e-02 1.874988e-03 3.327050e+00
    ##  [61] 4.376032e+00 3.381590e-02 5.139415e-01 1.410192e-02 1.368207e+00
    ##  [66] 3.121774e+00 3.137873e-01 3.435457e-01 1.166222e-01 3.573400e-01
    ##  [71] 2.474140e-01 4.460651e+00 1.334100e+00 2.410610e+00 6.611826e-02
    ##  [76] 1.710079e-01 3.349689e-01 2.619254e-01 2.444059e+00 1.388849e+00
    ##  [81] 4.011368e-01 1.543466e-01 1.641636e-02 3.093996e-02 1.067187e+00
    ##  [86] 1.074935e-01 2.028620e+00 5.629416e-03 4.565196e-01 3.765260e-02
    ##  [91] 4.384958e-01 1.172955e+00 8.148260e-06 3.541850e+00 1.690008e+00
    ##  [96] 8.246132e-01 7.644235e-01 1.812296e+00 2.033065e-02 7.437040e-01

#### fitted

This function should give you the fitted response variables. This is *not* the response variables in the data you fitted the model to, but instead the predictions that the model makes.

``` r
fitted.blm <- function(object, ...) predict(object, ...)

fitted(fit)
```

    ## $mean
    ##   [1]  0.1878164  2.3775693  2.5223758  4.0275633  5.0617376  6.4461740
    ##   [7]  6.4709281  7.9135870  8.8965125  9.5736535 11.0662204 11.8384177
    ##  [13] 13.0658631 13.9951290 15.0867391 15.9722627 17.0448547 17.9692453
    ##  [19] 18.9553066 20.0763959 21.0021000 21.9883158 22.9697553 23.8496716
    ##  [25] 24.9475380 25.9739151 26.8926936 27.9825421 28.8742717 29.9207586
    ##  [31] 31.0340463 31.8812860 32.9553394 33.8754448 34.9995193 35.9399402
    ##  [37] 36.9526243 37.9024226 38.9479576 39.9397711 40.9846418 41.9206737
    ##  [43] 42.9455939 43.9742365 44.9771352 45.9067469 46.9265292 47.9836278
    ##  [49] 48.9368328 49.8663580 50.9187970 51.8952884 52.9041011 53.9642850
    ##  [55] 54.8935605 55.9230430 56.9300483 57.9262931 58.9181972 59.8673100
    ##  [61] 60.8595170 61.9210760 62.8965444 63.9172011 64.9417967 65.9543247
    ##  [67] 66.9245067 67.8959477 68.9009407 69.9216558 70.9181683 71.8569199
    ##  [73] 72.9300019 73.8691123 74.9082706 75.9104492 76.8884298 77.8889564
    ##  [79] 78.8659633 79.8733950 80.8837781 81.9031450 82.8919191 83.8899864
    ##  [85] 84.9116731 85.8851160 86.9161775 87.8878236 88.8758315 89.8905757
    ##  [91] 90.8976954 91.9037442 92.8840310 93.9144961 94.8595653 95.8959168
    ##  [97] 96.8941673 97.8563591 98.8754006 99.8904133
    ## 
    ## $sd
    ##   [1] 1.136375 1.057428 1.069563 1.034658 1.034071 1.062177 1.072146
    ##   [8] 1.031417 1.030860 1.054341 1.028469 1.030061 1.026756 1.025009
    ##  [15] 1.025712 1.023383 1.023176 1.021846 1.021176 1.021833 1.019798
    ##  [22] 1.019061 1.018394 1.019969 1.017256 1.016614 1.016884 1.015580
    ##  [29] 1.016237 1.014823 1.014883 1.014592 1.013269 1.013870 1.012804
    ##  [36] 1.012215 1.011870 1.011933 1.011311 1.011077 1.011052 1.010739
    ##  [43] 1.010472 1.010466 1.010381 1.010250 1.010027 1.010268 1.009913
    ##  [50] 1.010580 1.009941 1.010139 1.010104 1.010277 1.010357 1.010321
    ##  [57] 1.010477 1.010648 1.010846 1.011473 1.011822 1.011581 1.011923
    ##  [64] 1.012185 1.012647 1.013158 1.013299 1.013698 1.014111 1.014599
    ##  [71] 1.015067 1.015846 1.016186 1.016748 1.017196 1.017803 1.018400
    ##  [78] 1.019037 1.019801 1.020430 1.021100 1.021862 1.022590 1.023369
    ##  [85] 1.024284 1.024998 1.026030 1.026731 1.027619 1.028566 1.029556
    ##  [92] 1.030581 1.031458 1.032732 1.033503 1.034663 1.035744 1.036750
    ##  [99] 1.037888 1.039138

#### plot

This function plots your model. You are pretty free to decide how you want to plot it, but I could imagine that it would be useful to see an x-y plot with a line going through it for the fit. If there are more than one predictor variable, though, I am not sure what would be a good way to visualise the fitted model. There are no explicit rules for what the `plot` function should do, except for plotting something so you can use your imagination.

``` r
# YOUR IMPLEMENTATION
```

#### predict

This function should make predictions based on the fitted model. Its interface is

``` r
predict(object, ...)
```

but the convention is that you give it new data in a variable `newdata`. If you do not provide new data, it instead gives you the predictions on the data used to fit the model.

``` r
predict.blm <- function(object, ...){

    mxy <- object$posterior$mu
    Sxy <- object$posterior$Sigma

    formula <- object$terms
    responseless.formula <- delete.response(terms(formula))
    frame <- model.frame(responseless.formula, ...)
    theta_x <- model.matrix(responseless.formula, frame)

    means = vector(length=nrow(theta_x))
    for(i in seq_along(means)) means[i] <- t(mxy) %*% theta_x[i,]
    sds = vector(length=nrow(theta_x))
    for(i in seq_along(sds)) sds[i] <- 1/object$beta + (t(theta_x[i,]) %*% Sxy %*% theta_x[i,])
    list(mean=means, sd=sds)
}

newdata <- data.frame(x=rnorm(5), z=rnorm(5))

predict(fit, data=newdata)
```

    ## $mean
    ## [1] -0.09712757 -0.11407059  4.43653088  3.13808140  2.43936659
    ## 
    ## $sd
    ## [1] 1.173537 1.043893 1.941064 3.093617 1.847414

#### print

This function is what gets called if you explicitly print an object or if you just write an expression that evaluates to an object of the class in the R terminal. Typically it prints a very short description of the object.

For fitted objects, it customarily prints how the fitting function was called and perhaps what the fitted coefficients were or how good the fit was. You can check out how `lm` objects are printed to see an example.

If you want to print how the fitting function was called you need to get that from when you fit the object in the `blm` constructor. It is how the constructor was called that is of interest, after all. Inside that function, you can get the way it was called by using the function `sys.call`.

``` r
print.blm <- function(x){
    cat("blm model: "); print(x$terms)
    cat("\n")
    cat("Posterior:\n")
    print(t(x$posterior$mu)); cat("\n")
    print(x$posterior$Sigma)
}

fit
```

    ## blm model: y ~ x + z
    ## 
    ## Posterior:
    ##      (Intercept)        x         z
    ## [1,]    1.561226 0.998953 -1.577903
    ## 
    ##               (Intercept)             x             z
    ## (Intercept)  0.3960735566 -4.926452e-04 -0.3622802757
    ## x           -0.0004926452  1.167065e-05 -0.0000910735
    ## z           -0.3622802757 -9.107350e-05  0.3675847517

#### residuals

This function returns the residuals of the fit. That is the difference between predicted values and observed values for the response variable.

``` r
# YOUR IMPLEMENTATION
```

#### summary

This function is usually used as a longer version of print. It gives you more information about the fitted model.

It does more than this, however. It returns an object with summary information. What that actually means is up to the model implementation so do what you like here.

``` r
# YOUR IMPLEMENTATION
```
