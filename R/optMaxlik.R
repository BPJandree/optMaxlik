devtools::use_package("maxLik")

#' Get number of unrestricted parameters from a maxLik object.
#'
#' Get number of unrestricted parameters from a maxLik object.
#' @param x maxLik object,
#' @keywords free parameters
#' @return numeric
#' @export


freePars <- function (x){
	pars <- coef(x)
	return(length(pars[pars!=0]))
}

#' Corrected AIC
#'
#' Compute the corrected AIC from a maxLik object.
#' @param x maxLik object
#' @param n number of observations for correction. Defaults to Inf, which coincides with the AIC.
#' @param penalty numeric that should be added to the log likelihood to adjust when using penalized likelihood or regularized estimators.
#' @keywords AICc
#' @return numeric
#' @export


aicc <- function (x,n=Inf,penalty=0){
	k=freePars(x) 
	c<-(2*k*(k+1))/(n-k-1)
	aic=2*k-2*(logLik(x)+penalty)
	aicc= aic + c 
	return(aicc)
}

#' Get the breads to build standard errors.
#'
#' Computes the inverse Fischer Information matrix
#' @param model maxLik object
#' @keywords Inverse Information matrix
#' @return matrix
#' @export

breads <-function(model){
  FischerI = -model$hessian 
  inv_I = solve(FischerI)
  return(inv_I)
}

#' Get standard errors.
#'
#' Computes the standard errors from the Hessian.
#' @param model maxLik object
#' @keywords standard errors
#' @return matrix
#' @export


se <- function (model){
  return(sqrt(diag(breads(model))))
}

#' Drop NA values from a vector.
#'
#' Removes NA values from a vector.
#' @param x vector
#' @keywords drop missing values
#' @return vector
#' @export


dropna <- function (x) {x[complete.cases(x)]}

#' Mainly auxiliary. Returns a vector indicating which parameters are fixed at zero.
#'
#' Returns a vector indicating which parameters are fixed at zero.
#' @param x vector
#' @keywords fixed parameters
#' @return vector
#' @export


areFixed <- function(x){
	pars <- abs(x)
	pars[pars==0]<-1
	pars[pars!=1]<-NA
	return(dropna(1:length(pars)*pars))
}


#' Mainly auxiliary. Inserts the content of an object into the supplied vector.
#'
#' Function to insert values into a vector.
#' @param a vector in which to insert a value
#' @param pos numeric indicating the position where to insert a value.
#' @param ... objects to be inserted.
#' @keywords fixed parameters
#' @return vector
#' @export
#' @examples
#' a <- 1:10
#' insert.at(a, 5, 0)

insert.at <- function(a, pos, ...){
	dots <- list(...)
	if(pos>1){
		pos=pos-1
	    stopifnot(length(dots)==length(pos))
	    result <- vector("list",2*length(pos)+1)
	    result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
	    result[c(FALSE,TRUE)] <- dots
	    unlist(result)
	    } else if (pos==1) {
	    	robj <- t(data.frame(c(c(dots),a)))
	    	names(robj)[1]<-"insert"
	    	if(length(names(a))==0){names(a)<-as.character(1:length(a))}
	    	names(robj)[-1]<-names(a)
	    	return(c(robj))
	    } else {
	    	stop("position should be >0")
	    }
}


#' Mainly auxiliary. Function to convert final parameters back into the format of the starting parameter vectors.
#'
#' Function to convert final parameters back into the format of the starting parameter vectors. Sets parameters that are contstrained in the final results to zero.
#' @param results vector of parameters
#' @param start vector of the original starting values.
#' @keywords reformat parameter estimates.
#' @return vector
#' @export



parsToOriginal <- function(results, start){
	start = start*0
	start[names(results)]<-results
	return(start)
}

#' Mainly auxiliary. Function to convert final parameters back into the format of the starting parameter vectors.
#'
#' Function to convert final parameters back into the format of the starting parameter vectors. Sets parameters that are contstrained in the final results to supplied starting values.
#' @param results vector of parameters
#' @param start vector of the original starting values.
#' @keywords reformat parameter estimates.
#' @return vector
#' @export


parsToOriginal2 <- function(results, start){
	start[names(results)]<-results
	return(start)
}

#' Mainly auxiliary. Same as as.numeric() but for matrices.
#'
#' Function to return a matrix with numericals.
#' @param mat matrix
#' @keywords as numeric
#' @return matrix
#' @export


as.numeric.matrix <- function(mat){
  F <- function(x){as.numeric(x)}
  return(apply(mat,2,F))
}

#' Main function of the package. Automated general to specific approach for custom Likelihood functions.
#'
#' Numerical minimization of AICc by optimization of a likelihood function and automated constraining of parameters.
#' @param LL a log likelihood function.
#' @param start named vector of starting values. 
#' @param initialfix vector indicating the parameters that should be treated as constants with values supplied in the start vector.
#' @param nobs number of observations to be used in the corrction of the AIC(c), defaults to Inf.
#' @param method numerical algorithm. See maxLik package. Defaults to BFGS
#' @param penalized vector indicating which parameters are penalized in your likelihood function. Will be used to compute adjustments for the AIC(c).
#' @param pw the weight of you penalty. Supported penalties are abs(parameter value)^pw.
#' @param nodrop vector indicating the parameters that will not be dropped. For example, distribution parameters of the ML function.
#' @return maxLik object of final model.
#' @export
#' @examples
#' library(maxLik)
#' library(compiler) 
#' # Fill a matrix with some random data.
#' mydata<-matrix(rnorm(1000), ncol=10, nrow=100)
#' 
#' # Create you own Log likelihood function.
#' crossectionalARMA_44 <- function(pars, ret="LL", Y){
#' 	 data=Y
#' 	 Log.L <-numeric()[1:T]
#' 	 b0    <-pars[1]
#' 
#' 	 B.y   <-pars[2]
#' 	 B.y2  <-pars[3]
#' 	 B.y3  <-pars[4]
#' 	 B.y4  <-pars[5]
#' 
#' 	 B.e   <-pars[6]
#' 	 B.e2  <-pars[7]
#' 	 B.e3  <-pars[8]
#' 	 B.e4  <-pars[9]
#' 
#' 	 nu    <-pars[10]
#' 	 sigma <-pars[11]
#' 
#' 	 df = as.numeric.matrix(data)
#' 	 T=nrow(df)
#' 	 N=ncol(df)
#' 
#' 	 p=4
#' 
#' 	 e = matrix(0,T,N)
#' 	 e[1:p,] <- 0 ## Initialize with e1 = 0
#' 
#' 	 Log.L[1:p]<-0
#' 
#' 	 A =N*log((gamma((nu+1)/2))/(((pi*(nu-2))^0.5)*gamma(nu/2))) - 0.5*N*log(max(0,sigma)^2)
#' 	 A2=((nu+1)/2)
#' 	 A3=(max(0,sigma)^2*(nu-2))  
#' 
#' 	 for (t in (p+1):T) {
#' 	  y = df[t,]#as.numeric(matrix(c(df[t,])))
#' 	  ymin =  df[t-1,]#as.numeric(matrix(c(df[t-1,])))
#' 	  ymin2 =  df[t-2,]#as.numeric(matrix(c(df[t-2,])))
#' 	  ymin3 =  df[t-3,]#as.numeric(matrix(c(df[t-3,])))
#' 	  ymin4 =  df[t-4,]#as.numeric(matrix(c(df[t-4,])))
#' 
#' 	  emin=e[t-1,]
#' 	  emin2=e[t-2,]
#' 	  emin3=e[t-3,]
#' 	  emin4=e[t-4,]
#' 
#'    MA = B.e*emin + B.e2*emin2 + B.e3*emin3 + B.e4*emin4 
#' 	  e[t,] = y - b0 - B.y*ymin - B.y2*ymin2 - B.y3*ymin3 - B.y4*ymin4 - MA
#' 
#' 	  }
#' 
#' 
#' 	 Log.L <- A - A2*(rowSums(log(1+ (e)^2 / A3)))
#' 	 Log.L[1:p]<-0
#' 
#' 	 if(ret=="e"){return(e)}else if (ret =="LLvec") {return(Log.L)} else(sum(Log.L))
#' 
#' } 
#'
#' # Optionally compile your function 
#' 
#' crossectionalARMA_44<-cmpfun(crossectionalARMA_44)
#' 
#' # opt.Maxlik takes only functions that have only a parameter vector as input.
#' # fix the data argument.
#' 
#' LL <- function(pars){crossectionalARMA_44(pars, ret="LL", Y=mydata)}
#' 
#' # create a vector of names parameter starting values.
#' 
#' start =c(b0=0, b1=0, b2=0, b3=0,b4=0, ma1=0, ma2=0, ma3=0, ma4=0, nu=120, sigma=sd(mydata))
#' 
#' # t-estimation
#' results <- opt.maxLik (LL=LL, start=start,initialfix=c(), nobs=dim(mydata)[1]*dim(mydata)[2])
#' summary(results)
#'
#' estimates=coef(results)
#'
#' parsToOriginal(estimates, start)
#'
#' # approximate normal estiamtion by fixing the degrees of freedom
#' results2 <- opt.maxLik (LL=LL, start=start,initialfix=c(10), nobs=dim(mydata)[1]*dim(mydata)[2])
#' summary(results2)
#' 
#' estimates2=coef(results2)
#' 
#' parsToOriginal(estimates2, start)
#' parsToOriginal2(estimates2, start)







 opt.maxLik <- function (LL, start, initialfix=numeric(), nobs=Inf, method="BFGS", penalized=numeric(), pw=2, nodrop=numeric()){



	newPars <- function (new.start, fix){

		activepars =rep(1, length(new.start))

		activepars[fix]<-0

		active.par.values=new.start[-fix]

		return(active.par.values)

	}



	maxID <- function(x){

		df<-cbind(1:length(x),x)

		df[df[,2]==max(x),1]

	}



	fixNew <- function(x, nodrop){

		old.fit = x

		old.tvec = coef(old.fit)/se(old.fit)

		old.tvec[is.na(old.tvec)] <-0

		new.freecoefs <- abs(1/old.tvec)

		new.freecoefs[new.freecoefs==Inf]<-0

		new.freecoefs[is.na(new.freecoefs)]<-0

		new.freecoefs[nodrop]<-0
		return(maxID(new.freecoefs))

	}



	newStart <- function(fit,fix){

		new.start <- coef(fit) 

		new.start[fix]<-0

		return(new.start)

	}



	ifnot <- function (x){if(x){return(FALSE)}else{return(TRUE)}}

		# insert names if needed
		if(length(names(start))!=length(start)){names(start)<- paste("par", as.character(1:length(start)), sep="")}
	 	# insert names if nonunique names
		if(length(unique(names(start)))!=length(start)){names(start)<- paste("par", as.character(1:length(start)), sep="")}
		

		# store supplied parameters

		allpars = start

		# adjust LL if there are initially fixed parameters

		if(length(initialfix)>0){	

			initialfixvals <- allpars[initialfix]

			start = start[-initialfix]

			use.LL <- function(start){

				runpars <-parsToOriginal(start,allpars)

				runpars[names(initialfixvals)]<-initialfixvals

				LL(runpars)

			}

		} else {

			use.LL <- LL

		}



	# initialize pointers

	fixed=numeric()

	iter=length(fixed) + 1



	# initial fit

	first.fit <- maxLik::maxLik(use.LL , start=start, method=method)

	penalty=sum(abs(parsToOriginal(coef(first.fit),start)[penalized]))^pw

	first.aic <- aicc(first.fit,nobs, penalty)


	notdrop = names(coef(first.fit)[nodrop])
	first.fixed <- fixNew(first.fit, notdrop)

	fixed[iter] <- first.fixed 



	if(ifnot(is.null(names(allpars)))){

		drop=names(start)[first.fixed]

	} else {

		drop=as.character(first.fixed)

	}

	message(paste("dropping",drop))



	# set next drop par to zero

	new.start <- newStart(first.fit, first.fixed)



	# remove drop par

	new.pars <- newPars(new.start, fixed)



	# build new LL wit hfixed par

	new.LL <- function(new.pars){

		runpars <-new.pars

		locs =rev(fixed)

		for (loc in 1: length(locs)){

			runpars<-insert.at(runpars, locs[loc], 0)

		}

		if(length(initialfix)>0){	

			runpars[names(initialfixvals)]<-initialfixvals

		}

		use.LL(runpars)

	}



	new.fit <- maxLik(new.LL, start=new.pars, method=method)

	penalty=sum(abs(parsToOriginal(coef(first.fit),start)[penalized]))^pw

	new.aic <- aicc(new.fit,nobs,penalty)


	if(length(new.aic) + length(first.aic) !=2){
		return(first.fit)
		stop("error in model convergence, returning first model fit. Are distribution parameters supplied in nodrop ?")
	}

	if(new.aic <= first.aic){continue=TRUE}else{continue=FALSE}



	if (continue){

		while (continue){



			iter = iter +1

			old.aic = new.aic

			old.fit = new.fit

	

			new.fixed <- fixNew(old.fit, notdrop) 

			fixed[iter] <- new.fixed 



			if(ifnot(is.null(names(allpars)))){

				drop=names(new.pars)[new.fixed]

			} else {

				drop=as.character(new.fixed)

			}

			message(paste("dropping",drop))

	

			new.start <- newStart(old.fit, new.fixed)



			new.pars <- newPars(new.start, new.fixed)



			new.LL <- function(new.pars){

				runpars <-new.pars

				locs =rev(fixed)

				for (loc in 1: length(locs)){

					runpars<-insert.at(runpars, locs[loc], 0)

				}

				if(length(initialfix)>0){	

					runpars[names(initialfixvals)]<-initialfixvals

				}

				use.LL(runpars)

			}



			new.fit <- maxLik(new.LL, start=new.pars, method="BFGS")

			penalty=sum(abs(parsToOriginal(coef(first.fit),start)[penalized]))^pw

			new.aic <- aicc(new.fit,nobs,penalty)

	

			if(new.aic <= old.aic){continue=TRUE}else{continue=FALSE}
	

		}


			return(old.fit)

		} else {

			return(first.fit)

		}

}
