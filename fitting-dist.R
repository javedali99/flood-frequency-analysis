# Term Project: Flood Frequency Analysis

library(readxl)
library(sqldf)
library(fitdistrib)
library(devtools)
library(EnvStats)
library(nsRFA)
library(MASS)
library(FAdist)
library(Lmoments)
library(evd)
#library(smwrBase)
library(fitdistrplus)
library(lmomco)
library(e1071)
library(goftest)
library(extRemes)
# -------------------------------------------------------------------------------------------
# IMPORT & PROCESS DATA
# -------------------------------------------------------------------------------------------
# read in data
#data <- read_excel('data/USGS02489500.xls', sheet = 1)
preprocessing <- function(data) {
    # sort by date and convert to dataframe
    data_processed <- data.frame(data[order(data$DecYear),])
    colnames(data_processed)[2] <- "aps"
    
    # remove period prior to 1900 (not enough data)
    data_processed <- data_processed[data_processed$Year > 1899,]
    
    # drop minimum value in years with multiple values
    data_processed <- sqldf("select Year, max(aps) as aps from data_processed group by Year")
    
    #-----------------------------------------------------------------------
    # Water Resources method:
    
    ## high outlier checking 
    # K_values <- read.csv('data/k_values.csv', header = TRUE)
    # K <- K_values[which(length(data_processed$aps) == K_values$n), "K"]
    # lmu <- mean(log(data_processed$aps)); ls <- sd(log(data_processed$aps))
    # upr <- exp(lmu + K*ls)
    
    # # remove outliers
    # data_cleaned <- data_processed[data_processed$aps < upr,]
    # 
    # # low outlier checking 
    # lmu <- mean(log(data_cleaned$aps)); ls <- sd(log(data_cleaned$aps))
    # lwr <- exp(lmu - K*ls)
    # 
    # # remove outliers
    # data_cleaned <- data_processed[data_processed$aps > lwr,]
    
    #-----------------------------------------------------------------------
    # IQR method
    q <- quantile(data_processed$aps)
    iqr <- q[4] - q[2]
    upr <- q[4] + 3*iqr
    lwr <- q[2] - 3*iqr
    outliers <- data_processed[(data_processed$aps > upr) | (data_processed$aps < lwr),]
    data_cleaned <- data_processed[(data_processed$aps < upr) & (data_processed$aps > lwr),]
    #-----------------------------------------------------------------------
    # # Standard deviation method: expect 99.7% of data to fall within this range
    # upr <- mean(data_processed$aps) + 3*sd(data_processed$aps)
    # lwr <- mean(data_processed$aps) - 3*sd(data_processed$aps)
    # data_cleaned <- data_processed[(data_processed$aps < upr) & (data_processed$aps > lwr),]
    #-----------------------------------------------------------------------
    # rank observations in decreasing order
    data_ranked <- data_cleaned[order(data_cleaned$aps, decreasing = TRUE),]
    data_ranked["rank"] <- seq(1,dim(data_cleaned)[1])
    data_ranked <- data_ranked[order(data_ranked$Year),]
    #attach(data_ranked)
    
    # order by rank
    data_sorted <- data_ranked[order(data_ranked$rank),]
    
    list("aps" = data_ranked$aps, 
         "ranked" = data_ranked, 
         "sorted" = data_sorted, 
         "processed" = data_processed,
         "noutliers" = dim(outliers)[1])
}

# for plotting
xseq <- seq(1,100000,10)

# -------------------------------------------------------------------------------------------
# ESTIMATION
# -------------------------------------------------------------------------------------------
# estimation for (log)normal distribution
normal <- function(x, method, log = FALSE, rp=NULL, plotx=FALSE) {
  n = length(x)
  if (log) { dist = "lognormal" } else { dist = "normal" }
  if ((method == "mle") || (method == "mme")) {
    mu <- mean(x); sd <-  sqrt(mean((x - mean(x))^2))
    params <- fitdistr(aps, dist)$estimate
    sds <- fitdistr(aps, dist)$sd
    mu <- params[1]; sd <- params[2]; 
  } else if (method == "pwme") { 
    lmom <- Lmoments(x)
    mu <- lmom[1] ; sd <- lmom[2]*sqrt(pi)
    sds <- c(sd/sqrt(n), 0.5113*sd/sqrt(n))
    }
  if (log) { 
    xt <- qlnorm(1-(1/rp), mu, sd)
    ploty <- dlnorm(plotx, mu, sd)
    lik <- sum(dlnorm(exp(x), mu, sd, log = TRUE))
    ad <- ad.test(exp(x), null = "plnorm", mean = mu, sd=sd)$statistic
    } else { 
    xt <- qnorm(1-(1/rp), mu, sd)
    ploty <- dnorm(plotx, mu, sd)
    lik <- sum(dnorm(x, mu, sd, log = TRUE))
    ad <- ad.test(x, null = "pnorm", mean = mu, sd=sd)$statistic
    }
  list("par" = list("mu" = mu, "sigma" = sd),
       "sd" = list("mu" = sds[1], "sigma" = sds[2]),
       "xt" = xt, 
       "ploty" = ploty, 
       "likelihood" = lik, 
       "parno" = 2, 
       "ad" = ad)
}
# -------------------------------------------------------------------------------------------
# estimation for  exponential distribution
exponential <- function(x, rp=FALSE, plotx=FALSE) {
  lambda <- 1/mean(x); n <- length(x)
  xt <- log(rp)*(1/lambda)
  ad <- ad.test(x, null = "pexp", rate=lambda)$statistic
  list("par" = list("lambda" = lambda), 
       "sd" = list("lambda" = sqrt(lambda^2 / n)),
       "xt" = xt, 
       "ploty" = dexp(plotx, lambda), 
       "likelihood" = sum(dexp(x, lambda, log=TRUE)), 
       "parno" = 1, 
       "ad" = ad)
}
# -------------------------------------------------------------------------------------------
# estimation for  gamma distribution
gam <- function(x, method, rp = NULL, plotx=FALSE) {
  if ((method == "mle")) {
    fit <- fitdist(x, "gamma", "mle" , start=list(shape=0.5, scale=0.5))
    params <- fit$estimate; sds <- fit$sd
    shape = params[1]; scale = params[2]
  } else if ((method == "mme")) {
    params <- fitdist(x, "gamma", "mme", start=list(shape=0.5, scale=0.5))$estimate
    shape = params[1]; scale = 1/params[2]
    sds <- NA
  } else if (method == "pwme") { 
    params <- gamma_Lmom(x)$estimate
    shape = params[1]; scale = 1/params[2]
    sds <- NA
  }
  ad <- ad.test(x, null = "pgamma", shape=shape, scale=scale)$statistic
  kt <- kt_gamma(x, rp)
  xt <- shape*scale + kt*sqrt(shape*(scale)^2)
  list("par" = list("shape" = shape, "scale" = scale), 
       "sd" = list("shape" = sds[1], "scale" = sds[2]),
       "xt" = xt, 
       "ploty" = dgamma(plotx, shape=shape, scale=scale), 
       "likelihood" = sum(dgamma(x, shape=shape, scale=scale, log=TRUE)), 
       "parno" = 2, 
       "ad" = ad)
}

# frequency factor for the gamma/Pearson 3 distribution
kt_gamma <- function(x, rp) {
  z <- qnorm(1 - (1/rp))
  k <- skew(x)/6
  kt <- z + (z^2 - 1)*k + (1/3)*(z^3-6*z)*k^2 - (z^2 - 1)*k^3+z*k^4+(1/3)*k^5
}
# -------------------------------------------------------------------------------------------
# estimation for Pearson 3 distribution
nll <- function(x, p) {
  alpha <- p[1]; beta <- p[2]; epsilon <- p[3]; n <- length(x)
  if(alpha>0 & epsilon <= min(x)) {
    -sum(dgamma(x-epsilon, scale = alpha, shape = beta, log=TRUE))
  } else { 
    return(Inf)
  }
}

p3 <- function(x, method, rp=NULL, plotx=FALSE) {
  if (method == "mle") {
    fit <- optim(x = x, c(500,3,200), nll, hessian=TRUE)
    params <- fit$par; hess <- fit$hessian
    location <- params[3]; scale <- params[1]; shape <- params[2]
    sds <- sqrt(abs(diag(solve(-hess))))
  } else if (method == "mme") { 
    params <- pearson_mom(x)$estimate
    location <- params[1]; scale <- params[2]; shape <- params[3]
    sds <- NA
  } else if (method == "pwme") { 
    params <- pearson_Lmom(x)$estimate
    location <- params[1]; scale <- params[2]; shape <- params[3]
    sds <- NA
  }
  kt <- kt_gamma(x, rp)
  xt <- location + shape*scale + kt*sqrt(shape*(scale)^2)
  ploty <- PearsonDS::dpearsonIII(plotx, scale = scale, shape = shape, location = location)
  lik <- sum(PearsonDS::dpearsonIII(x, scale = scale, shape = shape, location = location, log=TRUE))
  ad <- ad.test(x-location, null = "pgamma", shape=shape, scale=scale )$statistic
  list("par" = list("location" = location, "scale" = scale, "shape" = shape), 
       "sd" = list("location" = sds[3], "scale" = sds[1], "shape" = sds[2]),
       "xt" = xt, 
       "ploty" = ploty, 
       "likelihood" = lik, 
       "parno" = 3, 
       "ad" = ad)
}

# -------------------------------------------------------------------------------------------
# estimation for Log-Pearson 3 distribution
lp3 <- function(x, method, plotx=FALSE, rp=NULL) {
  if (method == "mle") {
    fit <- optim(x = x, c(3,3,4), nll, hessian=TRUE, method="BFGS")
    params <- fit$par; hess <- fit$hessian
    location <- params[3]; scale <- params[1]; shape <- params[2]
    sds <- sqrt(abs(diag(solve(-hess))))
  } else if (method == "mme") { 
    params <- pearson_mom(x)$estimate
    location <- params[1]; scale <- params[2]; shape <- params[3]
    sds <- NA
  } else if (method == "pwme") { 
    params <- parpe3(lmoms(x), checklmom = TRUE)$para
    mu <- params[1]; std <- params[2]; skw <- params[3]
    shape <- (2/skw)^2; scale <- std/sqrt(shape); location <- mu - scale*shape
    sds <- NA
  }
  kt <- kt_gamma(x, rp)
  xt <- exp(location + shape*scale + kt*sqrt(shape*(scale)^2))
  ad <- ad.test(x - location, null = "pgamma", shape = shape, scale = scale )$statistic
  list("par" = list("location" = location, "scale" = scale, "shape" = shape), 
       "sd" = list("location" = sds[3], "scale" = sds[1], "shape" = sds[2]),
       "ad" = ad,
       "xt" = xt,
       "ploty" = dgamma(log(plotx) - location, shape = shape, scale = scale)/plotx,
       "likelihood" = sum(dgamma(x - location, shape = shape, scale = scale, log = TRUE)) - sum(x), 
       "parno" = 3)
}
# -------------------------------------------------------------------------------------------
# estimation for Gumbel distribution
gumb <- function(x, method, rp = NULL, plotx=FALSE) {
  fit1 <- eevd(x, method = method, ci=TRUE, ci.parameter = "location")
  fit2 <- eevd(x, method = method, ci=TRUE, ci.parameter = "scale")
  params <-  fit1$parameter; df <- fit1$interval$dof; n <- fit1$interval$dof
  ub1 <- fit1$interval$limits[2];  ub2 <- fit2$interval$limits[2]
  sd_loc <- (ub1-params[1])/qt(0.975,df); sd_scale <- (ub2-params[2])/qt(0.975,df)
  xt <- params[1] - params[2]*log(-log(1-1/rp))
  ad <- try(ad.test(x, null = "pgumbel", alpha=params[1], scale=params[2])$statistic, silent=TRUE)
  if("try-error" %in% class(ad)) { ad <- ad.test(x, null = "pgumbel", loc=params[1], scale=params[2])$statistic } 
  list("par" = list("location" = params[1], "scale" = params[2]), 
       "xt" = xt, 
       "ploty" = evd::dgumbel(plotx, loc = params[1], scale = params[2]), 
       "likelihood" = sum(evd::dgumbel(x, loc = params[1], scale = params[2], log=TRUE)), 
       "parno" = 2, 
       "sd" = list("location" = sd_loc, "scale" = sd_scale),
       "ad" = ad)
}
# -------------------------------------------------------------------------------------------
# estimation for Weibull distribution
weibull <- function(x, method, rp=FALSE, plotx=FALSE) {
  if (method == "mle") {
    fit <- fitdist(x, dweibull3, start = list(shape = skew(x), scale = sd(x), thres = min(x)-100))
    params <- fit$estimate; sd <- fit$sd
    shape <- params[1]; scale <- params[2]; threshold <- params[3]
    sds <- list("shape" = sd[1], "scale" = sd[2], "threshold" = sd[3])
    ad <- ad.test(x, null = "pweibull3", shape=shape, scale=scale, thres = threshold)$statistic
    xt <- threshold + scale*(log(rp))^(1/shape)
  } else if (method == "mme") {
    params <- eweibull(x, method)$para
    shape = params[1];scale=params[2]
    ad <- ad.test(x, null = "pweibull", shape=shape, scale=scale)$statistic
    xt <- scale*(log(rp))^(1/shape)
    return(list("par" = list("shape" = shape, "scale" = scale),
         "xt" = xt,
         "ploty" = dweibull(plotx, shape = shape, scale = scale),
         "likelihood" = sum(dweibull(x, shape = shape, scale = scale, log=TRUE)),
         "parno" = 2,
         "sd" = NA,
         "ad" = ad))
  } else if (method == "pwme") {
      params <- parwei(lmoms(x),checklmom=TRUE)$para
      shape = params[3]; scale = params[2]; threshold = -params[1]
      ad <- ad.test(x, null = "pweibull3", shape=shape, scale=scale, thres = threshold)$statistic
      xt <- threshold + scale*(log(rp))^(1/shape)
      sds <- NA
  }
  if ((method=="mle") || (method=="pwme")) { 
    list("par" = list("shape" = shape, "scale" = scale, "threshold" = threshold),
         "xt" = xt,
         "ploty" = dweibull3(plotx, shape = shape, thres = threshold, scale = scale),
         "likelihood" = sum(dweibull3(x, shape = shape, thres = threshold, scale = scale, log=TRUE)),
         "parno" = 3,
         "sd" = sds,
         "ad" = ad)
   }
}

