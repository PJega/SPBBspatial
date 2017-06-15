# spbbspatial
Saddlepoint-Based Bootstrap (SPBB) Inference for Spatial Dependence

The SPBB method is used to make accurate (second-order) inference on a parameter of a model when 
there is no analytical solution for the distribution of an estimator.

Generally, the estimators are solution of some estimating equations. 
That's a particular estimator can be considered as a root of an estimating equation.

We proposed the SPBB method to make use of the relationship between a quadratic estimating equation (QEE)
in the normal random variables and the corresponding root to construct a confidence interval for the parameter of our interest.

We consider the SPBB method is an indirect method to construct a confidence interval for the parameter because 
this method constructs the distribution of the quadratic estimating equation using the saddlepoint approximations 
instead of directly approximating the distribution of the estimator.

The saddlepoint approximation for the distribution of QEE is derived by 
inverting the closed form expression of the moment generating function of QEE.

"spbbspatial" implements the SPBB method to make an inference of spatial dependence parmeter in the spatial 
regression models such as simultaneous autoregressive models (SAR), conditional autoregressive model, 
and simultaneous moving average model. 
By constructing the confidence interval for the spatial parameter in the regression model, 
we decide on which model should be used to reduce the prediction error.
