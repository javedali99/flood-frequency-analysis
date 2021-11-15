source('fitting-dist.R')

# turn off scientific notation
options(scipen = 1000)

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

# plot comparing different distributions
plot_combined <- function(station) {
  data <- read_excel(paste0("data/", station, ".xls"), sheet = 1) 
  out <- preprocessing(data)
  data_processed <- out$processed
  data_sorted <- out$sorted
  aps <- out$aps
  return_periods <- (dim(data_sorted)[1]+1)/data_sorted$rank
  
  best <- unlist(all_best[which(stations == station)])
  plot(log(return_periods), data_sorted$aps, xlab = "Return period (years)", ylab = "Discharge (cfs)", xaxt="n", 
       main = station)
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
  legend("bottomright", col = c(1, 2,3,8,5,6,4,"orange",9), lty = c(NA, rep(1,8)), pch=c(1, rep(NA,8)), c("Data", "Normal","Log-normal", "Exponential", "Gamma", "Pearson III", "Log-Pearson III", "Gumbel", "Weibull"))
  
  # chance ticks from log-space 
  at.x <- log(c(1,5,10,25,50,75,100))
  axis(1, at=at.x, labels=round(exp(at.x), 0), las=1)
}
