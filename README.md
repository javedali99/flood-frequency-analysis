# Flood frequency analysis of annual maximum discharge time series

Term project on "Flood Frequency Analysis" for "Water Resources in Changing Environment" class at [University of Central Florida](https://www.ucf.edu/).  

An interactive web application was built using the `Shiny` package in `R` for showing the results of analysis: [Shiny App](https://javedali.shinyapps.io/flood-frequency-analysis/)

### Data source
- [U.S. Geological Survey's (USGS) National Water Information System](https://nwis.waterdata.usgs.gov/usa/nwis/peak).
- Peak discharge data for 8 gauge stations on Pearl river, Mississippi were selected for the flood frequency analysis. 

<!---  
<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/15319503/141701862-c332913a-7072-4b11-a07a-7b6d8edbc505.png">
  <br>
    <em>Figure 1: Map of gauge stations locations on Pearl river</em>
</p>
--->

<div align="center">

*Table: USGS gauge stations*

|Station name | USGS Code |
|:-------: | :--------:  |
|Jackson | "USGS02486000" |
|Edinburg | "USGS02482000" |
|Carthage | "USGS02482550" |
|Lena | "USGS02483500" |
|Rockport | "USGS02488000" |
|Monticello | "USGS02488500" |
|Columbia | "USGS02489000" |
|Bogalusa | "USGS02489500" | 
  
</div>

### Probability distributions

We used following probability distributions for modelling annual maxima streamflow time series:
- Normal distribution
- Lognormal distribution
- Gamma distribution
- Pearson type 3 distribution
- Log-Pearson type 3 distribution
- Gumbel distribution
- Weibull distribution
- Exponential distribution

### Estimation methods

We selected following methods for estimating parameters of the distributions for this study: 
- Maximum Likelihood Estimation (MLE)
- Method of Moments (MOM)
- Probability Weighted Moments (PWM)

### Goodness of fit tests

Goodness-of-fit tests are used to summarise the discrepancy between a statistical model and the observed data. They are useful for comparing the observed values with either the values fitted by a model of interest or theoretical quantiles of a known sampling distribution. We used following metrics for determining whether a fit is satisfactory
or not. 
- Root-Mean-Square Error (RMSE)
- Kolmogorov-Smirnov test (K-S) 
- Anderson-Darling test (A-D) 
- Akaike Information Criterion (AIC)
- Bayesian Information Criterion (BIC)
- L-moment ratio diagram
