# Flood frequency analysis of annual maximum discharge time series

Term project on "Flood Frequency Analysis" for the "Water Resources in Changing Environment" class at the [Department of Civil, Environmental and Construction Engineering](https://www.cece.ucf.edu/) of [University of Central Florida](https://www.ucf.edu/).  

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

Maximum likelihood, ordinary moments and probability weighted moments (or L-moments) are widely used methods in frequency analysis to estimate parameters of the distributions that are also selected for this study.

- Maximum Likelihood Estimation (MLE)
- Method of Moments (MOM)
- Probability Weighted Moments (PWM)
