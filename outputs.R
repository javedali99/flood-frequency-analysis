# Author: Javed Ali

source('fitting.R')

stations <- list(Jackson = "USGS02486000",
                Edinburg = "USGS02482000",
                Carthage = "USGS02482550",
                Lena = "USGS02483500",
                Rockport = "USGS02488000",
                Monticello = "USGS02488500",
                Columbia = "USGS02489000",
                Bogalusa = "USGS02489500")

# best estimation method for each distribution
Jackson_best <- c("pwme","pwme", "mle", "mle", "mle", "mle","pwme", "pwme")
Bogalusa_best <- c("mle","mle", "mle", "mle", "mme", "pwme","pwme", "pwme")
Columbia_best <- c("mle","pwme", "mle", "pwme", "pwme", "pwme","pwme", "pwme")
Monticello_best <- c("mle","pwme", "mle", "mle", "pwme", "mle","pwme", "pwme")
Edinburg_best <- c("mle","pwme", "mle", "mle", "pwme", "mle","pwme", "pwme")
Carthage_best <- c("pwme","pwme", "mle", "mle", "pwme", "pwme","mle", "pwme")
Rockport_best <- c("pwme","pwme", "mle", "pwme", "pwme", "pwme","pwme", "pwme")
Lena_best <- c("pwme","pwme", "mle", "mle", "mle", "mle","mle", "pwme")
all_best <- list(Jackson_best, Edinburg_best, Carthage_best, Lena_best, Rockport_best, Monticello_best, Columbia_best, Bogalusa_best)


fit <- function(distr, method) {
  if (distr == "norm") {
    estim <- normal(aps, method, log = FALSE, rp = return_periods, plotx = xseq)
  }
  if (distr == "lognorm") {
    estim <- normal(log(aps), method, log = TRUE, rp = return_periods, plotx = xseq)
  }
  if (distr == "expo") {
    estim <- exponential(aps, rp = return_periods, plotx = xseq)
  }
  if (distr == "gam") {
    estim <- gam(aps, method, rp = return_periods, plotx = xseq)
  }
  if (distr == "p3") {
    estim <- p3(aps, method, rp = return_periods, plotx = xseq)
  }
  if (distr == "lp3") {
    estim <- lp3(log(aps), method, rp = return_periods, plotx = xseq)
  }
  if (distr == "gum") {
    estim <- gumb(aps, method, rp = return_periods, plotx = xseq)
  }
  if (distr == "wei") {
    estim <- weibull(aps, method, rp = return_periods, plotx = xseq)
  }
  estim
}
nstations <- length(stations)
lskew <- numeric(nstations)
lkurt <- numeric(nstations)
xskew <- numeric(nstations)
nobs <- numeric(nstations)

# turn off scientific notation
options(scipen = 1000)

for (i in 1:nstations) {
    station <- names(stations)[i]
    data <- read_excel(paste0("data/", stations[i], ".xls"), sheet = 1) 
    lmoment_ratios <- lmoms(data$`Annual Peak Streamflow (cfs)`)$ratios
    lskew[i] <- lmoment_ratios[3]
    lkurt[i] <- lmoment_ratios[4]
    
    pdf(file = paste0("Figures/", station, ".pdf"), 
        width     = 8,
        height    = 6)
    
    #png(file = paste0("Figures/", station, ".png"), 
    #    width     = 8,
    #    height    = 6,
    #    units     = "in",
    #   res       = 300,
    #   )
    
    
    out <- preprocessing(data)
    data_processed <- out$processed
    data_sorted <- out$sorted
    aps <- out$aps
    nobs[i] <- length(aps)
    return_periods <- (dim(data_sorted)[1]+1)/data_sorted$rank
    
    best <- unlist(all_best[i])
    plot(log(return_periods), data_sorted$aps, xlab = "Return period (years)", ylab = "Discharge (cfs)", xaxt="n")
    estim <- fit("norm", best[1])
    lines(log(return_periods), estim$xt, col=2, lwd=2)
    estim <- fit("lognorm", best[2])
    lines(log(return_periods), estim$xt, col=3, lwd=2)
    estim <- fit("expo", best[3])
    lines(log(return_periods), estim$xt, col=8, lwd=2)
    estim <- fit("gam", best[4])
    lines(log(return_periods), estim$xt, col=5, lwd=2)
    estim <- fit("p3", best[5])
    lines(log(return_periods), estim$xt, col=6, lwd=2)
    estim <- fit("lp3", best[6])
    lines(log(return_periods), estim$xt, col=4, lwd=2)
    estim <- fit("gum", best[7])
    lines(log(return_periods), estim$xt, col="orange", lwd=2)
    estim <- fit("wei", best[8])
    lines(log(return_periods), estim$xt, col=9, lwd=2)
    legend("bottomright", col = c(1, 2,3,8,5,6,4,"orange",9), lty = c(NA, rep(1,8)), pch=c(1, rep(NA,8)), 
           c("Data", "Normal","Log-normal", "Exponential", "Gamma", "Pearson III", "Log-Pearson III", "Gumbel", "Weibull"))
    at.x <- log(c(1,5,10,25,50,75,100))
    axis(1, at=at.x, labels=round(exp(at.x), 0), las=1)
    
    dev.off()
}

# log-Pearson 3 L-moments
# source: Griffis & Stedinger (2007)
tau0 <- c(-0.2308, -0.1643, -0.0740, 0, 0.0774, 0.1701, 0.2366)
coef <- matrix(c(0.0602, -0.1673, 0.8010, 0.2897, 0.0908, -0.1267, 0.7636, 0.2562, 0.1166, -0.0439, 0.6247, 0.2939, 0.1220, 0.0238,
         0.6677, 0.1677, 0.1152, 0.0639, 0.7486, 0.0645, 0.1037, 0.0438, 0.9327, -0.0951, 0.0776, 0.0762, 0.9771, -0.1394),
       byrow=TRUE, nrow=7)

tau4 <- matrix(NA, nrow=7, ncol=124)
for(i in 1:7) {
  tau3 <- seq(tau0[i], 1,0.01)
  if (124-length(tau3) > 0) {
    tau4[i,] <- c(rep(NA, 124-length(tau3)), coef[i,] %*% matrix(c(rep(1, length(tau3)),tau3, tau3^2, tau3^3), byrow=TRUE, nrow=4))
  } else {
    tau4[i,] <- coef[i,] %*% matrix(c(rep(1, length(tau3)),tau3, tau3^2, tau3^3), byrow=TRUE, nrow=4)
  }
}

# L-moment ratios diagram
png("Figures/L-moment-ratio-diagram.png", 
    width     = 8,
    height    = 6,
    units     = "in",
    res       = 300)

plotlmrdia(lmrdia(), xlim=c(0,1), xlab = "L-skewness", ylab ="L-kurtosis", noglo = TRUE, nosla = TRUE, nouni = TRUE, noray = TRUE, nocau = TRUE, nogpa = TRUE, nogov = TRUE, noaep4 = TRUE)

points(lskew, lkurt)
points(weighted.mean(lskew, nobs), weighted.mean(lkurt, nobs), pch=19, col="aquamarine3")

polygon(c(seq(tau0[1], 1,0.01),rev(seq(tau0[1], 1,0.01))), # LP3 region
        c(apply(tau4, 2, min, na.rm=TRUE), rev(apply(tau4, 2, max, na.rm=TRUE))),
        col=rgb(0, 0, 0.75,0.25), border = NA)

legend("topleft",
      col=c(2, 4, 2, 6, rgb(0, 0, 0.75,0.25), 2, 2, 8, 1, "aquamarine3"), 
      lty=c(NA, 2, NA, 1, 1, NA, 2, 1, NA, NA), 
      lwd=c(NA, 1, NA, 1, 8, NA, 1, 2, NA), NA, 
      pch=c(15, NA, 16, NA, NA, 17, NA, NA, 1, 19), 
      c("Normal", "Log-normal", "Exponential", "Pearson III", "Log-Pearson III", "Gumbel", "Weibull", "L-moment limit", "Stations", "Weighted mean"))

