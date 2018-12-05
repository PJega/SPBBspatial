        #############################################
        #INPUTS: z-vector of errors from SAR model  #
        #        W-spatial weight matrix            #
        #          This diagonal is zero,           #
        #                                           #
        #                                           #
        #############################################
require(MASS)
require(spdep)
require(psych)##Trace function

##START: ConfidenceIntervalSAR() to estimate MLE of rho, sigma square, conf.int for rho
ConfidenceIntervalSMA<-function(z,n=null,nrow=null,ncol=null,W=NULL,mn.r=NULL,mx.r=NULL){
        #z=spatial error vector -can get it from iterative least squares method or just remove the trend
        if(is.null(n)){
                n=nrow*ncol
        }else{n=n}
        
        
        if (is.null(W)&is.null(n)){
                #if the weight matrix is not given then the rook type weight matrix will be used
                W=nb2mat(cell2nb(nrow,ncol,type="rook"),style="W")       
        }else if(is.null(mn.r)|is.null(mx.r)){
                stop('Specify the max or min of spatial autocorrelatio based on weight matrix')
        }else{
                W=W     
        }
          
        
        I<-diag(seq(1,1,length=n))
        
        #default range for rho is based on the Binary weight matrix
        if(is.null(mn.r)){
                minrho=-0.25      
        }else{minrho=mn.r}
        
        if(is.null(mx.r)){
                maxrho=0.25      
        }else{maxrho=mx.r}
        
        #QEE to find MLE of rho
        qud.est.eqn=function(t){
                B=I+t*W
                B.inv=solve(B)
                Bt=I+t*t(W)
                B.Bt=B%*%Bt
                B.Btinv=solve(B.Bt)
                Q=tr(B.inv%*%W)*B.Btinv- 0.5*n*(B.Btinv%*%(W%*%Bt+B%*%t(W))%*%B.Btinv)
                val=t(z)%*%Q%*%z
                return(val)
        }
        
        #range for rhohat (mle)
        ranget=c(minrho,maxrho) 
        
        #compute mle of rho
        rhohat=uniroot(qud.est.eqn,ranget,tol=0.001)$root
       
         #compute mle of sigma square
        sigma.square.mle=(t(z)%*%solve((I+rhohat*W)%*%(I+rhohat*t(W)))%*%z/n)[1,1] 
        
        ##start: CDF.rho.mle()
        CDF.rho.mle=function(r0,t,sigma.square.mle)
        {   
                #Covaraince matrix of z
                sigma=function(r0)
                {
                        val=sigma.square.mle*(I+r0*W)%*%(I+r0*t(W))
                        return(val)   
                }
                
                #Qt  matrix
                Qt=function(t){ 
                        B=I+t*W
                        B.inv=solve(B)
                        Bt=I+t*t(W)
                        B.Bt=B%*%Bt
                        B.Btinv=solve(B.Bt)
                        Q=tr(B.inv%*%W)*B.Btinv- 0.5*n*(B.Btinv%*%(W%*%Bt+B%*%t(W))%*%B.Btinv)
                        return(Q)
                }
                
                
                sigma.r0.val=sigma(r0)
                Qt.val=Qt(t)
                pro.sigma.Qt=sigma.r0.val%*%Qt.val
                e.val=sort(Re(eigen(sigma.r0.val%*%Qt.val)$values))
                
                #To compute MGF
                Omega=function(s,pro.sigma.Qt){
                        val=I-2*s*pro.sigma.Qt
                        return(val)
                }
                
               
                ###Cumulant generating function of qud.est.eqn
                kgf=function(Omega.val)
                {
                        #log(mgf(s,t,r0))
                        val=-log(det(Omega.val))/2
                        return(val)
                        
                }
                
                
                #Second derivative of cumulant function
                kgfd2=function(Omega.val,pro.sigma.Qt)
                {
                        A=solve(Omega.val)%*%pro.sigma.Qt
                        val=2*tr(A%*%A)
                        return(val)
                }
                
                
                kgf.op=function(s,pro.sigma.Qt)
                {
                        #log(mgf(s,t,r0))
                        val=-log(det(Omega(s,pro.sigma.Qt)))/2
                        return(val)
                        
                }
                #compute the s hat
                sh=suppressWarnings(optimize(kgf.op,
                         interval=c(1/(2*e.val[1])+0.00001,1/(2*e.val[n])-0.00001),
                         pro.sigma.Qt=pro.sigma.Qt,tol=.001)$minimum)
               
                 #Define Wt function
                Wt=function(s,kgf.val)
                {
                        val=sign(s)*sqrt(-2*kgf.val)
                        return(val)
                }
                
                #Define Ut function
                Ut=function(s,kgf.val)
                {
                        val=s*sqrt(kgfd2.val)
                        return(val)
                }
               
                
                Omega.val=Omega(sh,pro.sigma.Qt)
                kgf.val=kgf(Omega.val)
                kgfd2.val=kgfd2(Omega.val,pro.sigma.Qt)
                Wt.val=Wt(sh,kgf.val)
                Ut.val=Ut(sh,kgfd2.val)
                
                #saddle CDF of MLE
                1-(pnorm(Wt.val)+dnorm(Wt.val)*((1/Wt.val)-(1/Ut.val)))
                
                
        }##END: CDF.rho.mle()
        
        ##START: CDF.rho.mle.tval() to get the function of rho
        CDF.rho.mle.tval=function(r0)
        {
                val=CDF.rho.mle(r0,rhohat,sigma.square.mle)
                return(val)
        }##END: CDF.rho.mle.tval()
        
        
        #CDF is decreasing function in rho
        #Inverting (Pivoting) CDF to find the confidence interval for rho
        
        ##START: Find.rho.lower.fun()
        Find.rho.lower.fun=function(rhol)
        {
                val=abs(CDF.rho.mle.tval(rhol)-0.975)
                return(val)
        }##END: Find.rho.lower.fun()
        
        rho.lower=function(){
                val=suppressWarnings(optimize(Find.rho.lower.fun,c(ranget[1],rhohat),tol=.0001)$minimum)
                return(val)
        }
        
        #CDF is decreasing function in rho
        Find.rho.upper.fun=function(rhou)
        {
                val=abs(CDF.rho.mle.tval(rhou)-0.025)
                return(val)
        }
        
        rho.upper=function(){
                val=suppressWarnings(optimize(Find.rho.upper.fun,c(rhohat,ranget[2]),tol=.0001)$minimum)
                return(val)
                
        }
        r.lower=rho.lower()
        r.upper=rho.upper()
        list.out=list(rhohat=rhohat,sigma.sq=sigma.square.mle,CI=c(r.lower,r.upper))
        return(list.out)
       
}##END: ConfidenceIntervalSAR()       


