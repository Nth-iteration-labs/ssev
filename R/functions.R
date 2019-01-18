#' Compute sample size
#'
#' Function to compute the optimal sample size for a comparison of two  means
#' (with equal or unequal variances) or proportions. Function returns the standard
#' sample size for an RCT with the specified power, as well as the optimal sample size
#' for a population of size N.  
#'
#' @param means A vector of length 2 containing the (assumed) means of the two groups
#' @param sds A vector containing the (assumed) standard deviations of the two groups. 
#' When only one element is supplied equal variances are assumed.
#' @param proportions A vector of length 2 containing the (assumed) proportions of the two groups
#' @param N Estimated population size
#' @param power Desired power for the classical RCT
#' @param sig.level Significance level of the test used (alpha)
#' @param ties Probability of choosing the first group in case of a tie (i.e., H0 is not rejected)
#' @param .verbose Whether or not verbose output should be provided, default FALSE
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of type ssev
#' @export
#'
#' @examples compute_sample_size(means=c(0,1), sds=2, N=100)
#' @examples compute_sample_size(means=c(0,1), sds=2, N=10000, power=.9)
#' @examples compute_sample_size(means=c(0,1), sds=c(1,2), N=10000)
#' @examples compute_sample_size(proportions=c(.5,.7), N=5000)
compute_sample_size <- function(
	means = NULL,					
	sds = NULL,
	proportions = NULL,
	N = Inf,
	power = .8,
	sig.level = .05,
	ties = .5,
	.verbose = FALSE, ...
	){
		
	if(length(means) > 2){
		stop("More than k=2 means not implemented.")
	}
	if(length(proportions) > 2){
		stop("More than k=2 proportions not implemented.")
	}
	
	## Compare means
	if(length(means) == 2){	
		if(.verbose){print("Comparing means.")}
		
		# Checks:
		if(length(sds)==0){
			stop("Unable to compare means without specifying standard deviations.")
		}
		if(length(means)!=length(sds) && length(sds) != 1){
			stop("Length of means and standard devaiations does not match.")
		}		
	
		# Equal variances:
		if(length(sds)==1){
			if(.verbose){print("Equal variances assumed.")}
						
			delta <- means[1]-means[2]
						
			# Classical:
			n_rct <- ceiling(MESS::power_t_test(
				power=power, 
				sd=sds[1], 
				delta=delta, 
				sig.level=sig.level,
				alternative="two.sided",
				type="two.sample")$n[1]
				)
						
			ev_rct <- ev_means_equal(n_rct, N, means, sds[1], sig.level, ties)
						
			# Optimal:
			n_opt <- ceiling(
				stats::optimize(ev_means_equal, c(2, N/2), 
				maximum=TRUE, N=N, means=means, sd=sds[1], sig.level=sig.level, ties=ties)$maximum
				)
			ev_opt <- ev_means_equal(n_opt, N, means, sds[1], sig.level, ties)
						
			method <- "Two independent samples t-test (equal variances assumed)."
			note <- "n is number in *each* group."
					
		# Unequal variances:
		} else {
			if(.verbose){print("Unequal variances assumed.")}
						
			delta <- means[1]-means[2]
			sd.ratio <- sds[2]/sds[1]
						
			# Classical:
			n_rct <- ceiling(MESS::power_t_test(
				power=power, 
				sd=sds[1], 
				delta=delta, 
				sig.level=sig.level,
				sd.ratio=sd.ratio,
				alternative="two.sided",
				df.method="welch",
				type="two.sample")$n[1]
				)
						
			ev_rct <- ev_means_unequal(n_rct, N, means, sds, sig.level, ties)
						
			# Optimal:
			n_opt <- ceiling(
				stats::optimize(ev_means_unequal, c(2, N/2), 
				maximum=TRUE, N=N, means=means, sds=sds, sig.level=sig.level, ties=ties)$maximum
				)
			ev_opt <- ev_means_unequal(n_opt, N, means, sds, sig.level, ties)
						
			method <- "Two independent samples t-test (unequal variances)."
			note <- "n is number in *each* group; using Welch test."

		}
		
		if(n_rct > N/2){
			note <- paste(note, "RCT size larger than population.")
		}
		
		gain <- (100*ev_opt)/ev_rct - 100
		
		# Nicer Nrct max half:
		if(n_rct * 2 > N){
			n_rct <- N/2
		}
		
		obj <- list(
			"method" = method,
			"note" = note,
			"n_rct" = n_rct,
			"ev_rct" = ev_rct,	
			"n_opt" = n_opt,
			"ev_opt" = ev_opt,
			"gain" = gain,
		     "params" = c(as.list(environment()), list(...))
			)
		class(obj) <- c("ssev", class(obj))
		return(obj)
				
	}
					
	## Compare proportions
	if(length(proportions)>1){
		if(.verbose){print("Comparing proportions.")}
		
		if(length(proportions)!=2){
			stop("Proportions should be a vector of length 2.")
		}
					
		# Classical:
		n_rct <- ceiling(MESS::power_prop_test(p1 = proportions[1], p2 = proportions[2], sig.level = sig.level,
		  power = power, ratio = 1, alternative = "two.sided", tol = .Machine$double.eps^0.25)$n[1])
					
		ev_rct <- ev_proportions(n_rct, N, proportions, sig.level, ties)
					
		# Optimal:
		n_opt <- ceiling(
			stats::optimize(ev_proportions, c(n_rct*.1, N/2), 
			maximum=TRUE, N=N, proportions=proportions, sig.level=sig.level, ties=ties)$maximum
			)
		ev_opt <- ev_proportions(n_opt, N, proportions, sig.level, ties)
					
		method <- "Comparing two proportions."
		note <- "n is number in *each* group."

		if(n_rct > N/2){
			note <- paste(note, "RCT size larger than population.")
		}
		
		gain <- (100*ev_opt)/ev_rct - 100
		
		# Nicer Nrct max half:
		if(n_rct * 2 > N){
			n_rct <- N/2
		}
				
		obj <- list(
			"method" = method,
			"note" = note,
			"n_rct" = n_rct,
			"ev_rct" = ev_rct,	
			"n_opt" = n_opt,
			"ev_opt" = ev_opt,
			"gain" = gain,
		    "params" = c(as.list(environment()), list(...))
			)
		class(obj) <- c("ssev", class(obj))
		return(obj)
				
	}
			
	warning("No valid cases found: please supply the appropriate arguments.")
	return(FALSE)
}
	



	
#' Compute expected value as function of n, N
#'
#' Comparing means with equal variances
#'
#' @param n Sample size per group
#' @param N Population size (estimate)
#' @param means Vector of estimated means
#' @param sd Standard deviation of the groups (assumed equal)
#' @param sig.level Significance level
#' @param ties Tie-breaking probability
#'
#' @return A scalar indicating the expected mean reward per unit in the population
ev_means_equal <- function(n, N, means, sd, sig.level, ties){
	
	if(N==Inf){
		stop("Not possible to compute population expected value for infinite populations.")
	}
	
	sample_size <- min(2*n,N)
	ev_experiment <- (.5*means[1] + .5*means[2])*sample_size
	
	sig <- sig.level/2
	delta <- means[1]-means[2]
	
	#if equal:
	pM1 <- MESS::power_t_test(n=n, sd=sd, delta=delta, 
				sig.level=sig,alternative="one.sided",type="two.sample")$power
	pM2 <- MESS::power_t_test(n=n, sd=sd, delta=-delta, 
				sig.level=sig,alternative="one.sided",type="two.sample")$power
		
	pH0 <- 1 - (pM1+pM2)
	
	ev_guideline <- (N - sample_size) * (
			(ties*pH0 + pM1) * means[1] +
			((1-ties)*pH0 + pM2) * means[2])

	return( (ev_experiment+ev_guideline)/N )
	
}




#' Compute expected value as function of n, N
#'
#' Comparing means with unequal variances
#'
#' @param n Sample size per group
#' @param N Population size (estimate)
#' @param means Vector of estimated means
#' @param sds Vector of standard deviation of the groups
#' @param sig.level Significance level
#' @param ties Tie-breaking probability
#'
#' @return A scalar indicating the expected mean reward per unit in the population
ev_means_unequal <- function(n, N, means, sds, sig.level, ties){

	if(N==Inf){
		stop("Not possible to compute population expected value for infinite populations.")
	}
	
	sample_size <- min(2*n,N)
	ev_experiment <- (.5*means[1] + .5*means[2])*sample_size
	
	sig <- sig.level/2
	delta <- means[1]-means[2]
	sd.ratio <- sds[2]/sds[1]
	sd <- sds[1]
	
	pM1 <- MESS::power_t_test(n=n, sd=sd, power=NULL, ratio=1, sd.ratio=sd.ratio, delta=delta, 
		sig.level=.025,alternative="one.sided",type="two.sample", df.method="welch")$power
	pM2 <- MESS::power_t_test(n=n, sd=sd, power=NULL, ratio=1, sd.ratio=sd.ratio, delta=-delta, 
			sig.level=.025,alternative="one.sided",type="two.sample", df.method="welch")$power
		
	pH0 <- 1 - (pM1+pM2)
	
	ev_guideline <- (N - sample_size) * (
			(ties*pH0 + pM1) * means[1] +
			((1-ties)*pH0 + pM2) * means[2])

	return( (ev_experiment+ev_guideline)/N )
	
}



#' Compute expected value as function of n, N
#'
#' Comparing proportions 
#'
#' @param n Sample size per group
#' @param N Population size (estimate)
#' @param proportions Vector of two proportions
#' @param sig.level Significance level
#' @param ties Tie-breaking probability
#'
#' @return A scalar indicating the expected mean reward per unit in the population
ev_proportions <- function(n, N, proportions, sig.level, ties){
	
	if(N==Inf){
		stop("Not possible to compute population expected value for infinite populations.")
	}
	
	sig <- sig.level/2
	sample_size <- min(2*n,N)
	ev_experiment <- (.5*proportions[1] + .5*proportions[2])*sample_size
	
	h <- 2*asin(sqrt(proportions[1]))-2*asin(sqrt(proportions[2]))
	
	pM1 <- tryCatch({
	    pwr::pwr.2p.test(h = h, n = n, sig.level = .05, alternative="greater")$power
	}, warning = function(w) {
	   	warning(paste('Power calculation unable to converge for n=',n,sep=""))
		return(0)
	})
	
	pM2 <- tryCatch({
	    pwr::pwr.2p.test(h = h, n = n, sig.level = .05, alternative="less")$power
	}, warning = function(w) {
	   	warning(paste('Power calculation unable to converge for n=',n,sep=""))
		return(0)
	})
	
	pH0 <- 1 - (pM1+pM2)
	
	ev_guideline <- (N - sample_size) * (
			(ties*pH0 + pM1) * proportions[1] +
			((1-ties)*pH0 + pM2) * proportions[2])

	return( (ev_experiment+ev_guideline)/N )	
}


#' Pretty printing of ssev object
#'
#'
#' @param x Object of type ssev for pretty printing
#' @param digits Standard number of digits for pretty printing, default is getOption("digits")
#' @param ... further arguments passed to or from other methods.
#'
#' @return NULL
print.ssev <- function(x, digits = getOption("digits"), ...){
    cat("\n")
	cat(x$method, "\n\n")
    note <- x$note
    x[c("method", "note")] <- NULL
    
	cat(paste(format("Sample size RCT", width = 35L, justify = "right"), 
        format(x$n_rct, digits = digits), sep = " = "), sep = "\n")
	cat(paste(format("Expected mean reward RCT", width = 35L, justify = "right"), 
	    format(x$ev_rct, digits = digits), sep = " = "), sep = "\n")
	cat(paste(format("Sample size optimal", width = 35L, justify = "right"), 
		format(x$n_opt, digits = digits), sep = " = "), sep = "\n")
	cat(paste(format("Expected mean reward optimal", width = 35L, justify = "right"), 
		format(x$ev_opt, digits = digits), sep = " = "), sep = "\n")
    
	cat("\n")
	cat(paste(format("Percentage gain (optimal over RCT):", width = 35L, justify = "right"), 
		format(x$gain, digits = digits), sep = " = "), sep = "\n")
	
	cat("\n")
	cat(paste("Sig level: ",x$params$sig.level,", power (RCT): ",x$params$power, ", population size (optimal): ", x$params$N, sep=""))

    cat("\n", "NOTE: ", note, "\n\n", sep = "")
    invisible(x)
}