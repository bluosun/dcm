#' @title Discrete Choice Model (DCM) for Nonignorable Missing Data
#' 
#' @description
#' Implements inverse probability weighting (IPW) with the logit discrete choice nonresponse model (LDCM).
#'
#' @details
#' This function implements IPW using LDCM weights as described in Tchetgen Tchetgen et al. (2017). The implementation is based on 
#' a default linear main effects model of the observed variables in the \emph{r}-th missing data pattern for the log ratio 
#' of probabilities for observing the \emph{r}-th missing data pattern versus the complete data.
#'
#'
#' @references
#' Tchetgen Tchetgen, E., Wang, L. and Sun, B. (2017). \href{http://www3.stat.sinica.edu.tw/ss_newpaper/SS-2016-0325_na.pdf}{Discrete Choice Models for Nonmonotone
#' Nonignorable Missing Data: Identification and Inference}. Statistica Sinica (doi: 10.5705/ss.202016.0325).
#' 
#' @param data_frame A data frame containing the variables in the regression model.
#' @param regr_formula An object of class "formula" describing the regression model to be fitted
#' @param regr_fam A description of the error distribution and link function to be used in the regression model. 
#' This can be a character string (e.g. "binomial" or "gaussian") naming a family function, same as the input to glm().
#'
#' @return A "dcm" object containing the following items:
#' \item{CC}{An object of class "glm" from complete-case regression analysis.}
#' \item{IPW}{The corresponding coefficient point estimate, bootstrap standard error, confidence interval and p-value from IPW analysis.}
#' \item{DAT}{A data frame containing the variables in the regression model and the missing data pattern indicator R.}
#'
#' @examples
#' 
#' #The 'airquality' dataset is included in R base and contains nonmonotone missing values
#' #in the variables Ozone (ppb) and Solar.R (lang). 
#' out <- dcm(airquality, Temp~Ozone+Solar.R+Wind, gaussian)
#'
#'@export




dcm <- function(data_frame, regr_formula, regr_fam) {

	if (is.matrix(data_frame)) {
		data_frame = data.frame(data_frame);
	}
	regr_frame = get_all_vars(regr_formula, data_frame);

	#following adapted from Linbo Wang's code
      ### 0. Descriptive ###
      #   percentage of missing
      N = nrow(regr_frame);
      pct.miss = apply(apply(regr_frame,2,is.na),2,sum)*100/N;
      ### 1. Complete Case analysis ###
      CC.result = glm(as.formula(regr_formula),family=regr_fam, data=regr_frame);
      ### 2. DCM: IPW
      #   mising data pattern R
      R2 = apply(regr_frame,1,function(x) sum(is.na(x)*(2^(0:(ncol(regr_frame)-1))))) + 1
	print(paste("%complete-cases: ",round(sum(R2==1)/N*100,2),sep=""))
	
	R=numeric(N); j=0;
	for (i in sort(unique(R2))) {j=j+1; R[R2==i]=j}

      #   complete case weights pi1.com
	ipw.est = ipw_dcm(regr_frame, R, regr_formula, regr_fam);
	print(ipw.est)
      #   bootstrap estimate of variance
      btstraps = matrix(NA,100,length(ipw.est))
      for(btstrap in 1:100){
            index = sample(1:N,N,replace=TRUE)
            R.b = R[index]; regr_frame.b = regr_frame[index,];
            btstraps[btstrap,] = ipw_dcm(regr_frame.b, R.b, regr_formula, regr_fam);
       }
	ipw.sd = sqrt(apply(btstraps, 2, var, na.rm=TRUE));
	ipw.z = ipw.est/ipw.sd; 
	alpha=0.05;
	ci_lower = ipw.est -stats::qnorm(1-alpha/2)*ipw.sd;
	ci_upper = ipw.est +stats::qnorm(1-alpha/2)*ipw.sd;
      ipw.result = rbind(ipw.est,ipw.sd,ci_lower,ci_upper,2*pnorm(-abs(ipw.z)));
      row.names(ipw.result)= c("Estimate","Std. Error","95% CI lower","95% CI upper", "p-val" );
	object <- list(CC=CC.result, IPW=t(ipw.result), DAT=cbind(R,regr_frame));
      class(object) <- "dcm";
	return(object);
}

ipw_dcm <- function(data_frame, R, regr_formula, regr_fam) {

      data = as.matrix(cbind(rep(1,nrow(data_frame)),data_frame))
      odds.com = matrix(NA,sum(R==1),max(R))
      odds.com[,1] = rep(1,sum(R==1))
      for(r in 2:max(R)){
          Rr = as.numeric(R==r);  R1r = R==1 | R==r
          if(sum(Rr)>1){
              obs.index = which(!is.na(data[R==r,][1,]))   
          }else if (sum(Rr)==1){
              obs.index = which(!is.na(data[R==r,][1]))  
          }else {
		  break;
	    }
          alphar = glm(Rr ~ data[,obs.index] - 1, subset = R1r, family=binomial)$coef
          odds.com[,r] = exp(as.matrix(data[R==1,obs.index]) %*% alphar) 
      }
      pi1 = 1/(rowSums(odds.com, na.rm=TRUE));
      environment(regr_formula) <- environment();
      ipw.result = suppressWarnings(glm(as.formula(regr_formula),family=regr_fam, weights=1/pi1,
                       data=data_frame[R==1,]));
	return(ipw.result$coef);
}

