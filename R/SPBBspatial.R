#' Saddlepoint-Based Bootstrap (SPBB) Inference for Spatial Dependence Parameter
#'
#' Saddlepoint-Based Bootstrap (SPBB) Inference for Spatial Dependence Parameter in
#' the Spatial Regression Models such as Simultaneous Autoregressive (SAR) Model,
#' Conditional Autoregressive (CAR) Model, and Simultaneous Moving Average (SMA) Model.
#'
#' @param formula regression formula in R
#' @param data list containing
#' @param nbw see \code{poly2nb}
#' @param family either 'SAR', or 'CAR', or 'SMA'
#' @param mn.r minimum value of spatial dependence parameter for the SMA model
#' @param mx.r maximum value of spatial dependence parameter for the SMA model
#' @return A list of coefficients of the spatial regression covariates and the confidence interval for the spatial dependence parameter
#' @author Pratheepa Jeganathan
#' @details This function takes the regression formula, weight matrix, and the name of the spatial
#' regression model (SAR, CAR, SMA) and calculates the spatial correlated residuals from the \code{spautolm}.
#' Then, this function returns the confidence interval for the spatial dependence
#' parameter using the SPBB method.
#' @seealso \code{spautolm}
#' @export
#' @source SAR_SPBB.R
#' @source CAR_SPBB.R
#' @source SMA_SPBB.R


library(spdep)
source('./R/SAR_SPBB.R')
source('./R/CAR_SPBB.R')
source('./R/SMA_SPBB.R')
SPBBspatial<-function(formula, data = list(), nbw, family='SAR',mn.r=-.9,mx.r=.9){
    #nbw:readPolyshap->poly2nb() should be nbw

        if(family=='SAR'){
                listw=nb2listw(nbw, style="W")
                sarfit=spautolm(formula=formula, data=data, listw=listw,family=family)
                betahat=matrix(sarfit$fit$coefficients,nrow= length(sarfit$fit$coefficients),ncol=1)
                mt <- terms(formula, data = data)
                mf <- lm(formula, data,method = "model.frame")
                y <- model.response(mf, "numeric")
                if (any(is.na(y))) stop("NAs in dependent variable")
                X<- model.matrix(mt, mf)
                if (any(is.na(X)))
                    stop("NAs in independent variable")
                if (NROW(X) != length(listw$neighbours))
                    stop("Input data and neighbourhood list have different dimensions")

                zhat=y-X%*%betahat
                W=nb2mat(nbw,style="W")
                z=as.vector(zhat)
                res=list(fit.s=summary(sarfit),SPBB=ConfidenceIntervalSAR(zhat,n=length(y),W=W))
                return(res)
        }else if(family=='CAR'){
            listw=nb2listw(nbw, style="B")
            sarfit=spautolm(formula=formula, data=data, listw=listw,family=family)
            betahat=matrix(sarfit$fit$coefficients,nrow= length(sarfit$fit$coefficients),ncol=1)
            mt <- terms(formula, data = data)
            mf <- lm(formula, data,method = "model.frame")
            y <- model.response(mf, "numeric")
            if (any(is.na(y))) stop("NAs in dependent variable")
            X<- model.matrix(mt, mf)
            if (any(is.na(X)))
                stop("NAs in independent variable")
            if (NROW(X) != length(listw$neighbours))
                stop("Input data and neighbourhood list have different dimensions")

            zhat=y-X%*%betahat
            W=nb2mat(nbw,style="B")
            z=as.vector(zhat)
            res=list(fit.s=summary(sarfit),SPBB=ConfidenceIntervalCAR(zhat,n=length(y),W=W))
            return(res)
        }else if(family=='SMA'){
            listw=nb2listw(nbw, style="W")
            sarfit=spautolm(formula=formula, data=data, listw=listw,family=family)
            betahat=matrix(sarfit$fit$coefficients,nrow= length(sarfit$fit$coefficients),ncol=1)
            mt <- terms(formula, data = data)
            mf <- lm(formula, data,method = "model.frame")
            y <- model.response(mf, "numeric")
            if (any(is.na(y))) stop("NAs in dependent variable")
            X<- model.matrix(mt, mf)
            if (any(is.na(X)))
                stop("NAs in independent variable")
            if (NROW(X) != length(listw$neighbours))
                stop("Input data and neighbourhood list have different dimensions")

            zhat=y-X%*%betahat
            W=nb2mat(nbw,style="W")
            z=as.vector(zhat)
            res=list(fit.s=summary(sarfit),SPBB=ConfidenceIntervalSMA(zhat,n=length(y),W=W,mn.r=mn.r,mx.r=mx.r))
            return(res)
        }
}
#source this file
#SPBBspatial(A ~ towns + pale, data=eire, eire.nb,family='SAR')
