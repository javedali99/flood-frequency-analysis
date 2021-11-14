library(shiny)
source('fitting.R')
source('combined-plots.R')

# Define UI for application that draws a histogram
ui <- navbarPage("Term Project: Flood Frequency Analysis",
     tabPanel("Explore Data",
     
       sidebarLayout(
          sidebarPanel(
            # choose gauging station
            selectInput("station", "Station",
                        c("Jackson" = "USGS02486000", 
                          "Edinburg" = "USGS02482000", 
                          "Carthage" = "USGS02482550", 
                          "Lena" = "USGS02483500", 
                          "Rockport" = "USGS02488000", 
                          "Monticello" = "USGS02488500", 
                          "Columbia" = "USGS02489000",
                          "Bogalusa" = "USGS02489500")),
            
            # select plot type
            selectInput("plottype", "Plot type",
                        c("Time series", 
                          "Density")),
              h5("Number of outliers detected:"),
              uiOutput('outliers')
          ),
          
          mainPanel(
             plotOutput("timePlot")
          )
        )
     ),
     tabPanel("Fit Distributions",
          sidebarLayout(
            sidebarPanel(
              # choose distribution
              selectInput("distr", "Distribution",
                          c("Normal" = "norm", 
                            "Log-normal" = "lognorm",
                            "Exponential" = "expo",
                            "Gamma" = "gam",
                            "Pearson 3" = "p3",
                            "Log-Pearson 3" = "lp3",
                            "Gumbel" = "gum",
                            "Weibull" = "wei")),
              
              # choose estimation method
              radioButtons("method",
                      "Estimation method",
                      c("Maximum Likelihood (MLE)" = "mle", 
                        "Method of Moments (MOM)" = "mme", 
                        "Probablity weighted moments (PWM)" = "pwme")),
              h5("RMSE:"),
              uiOutput('rmse'),
              h5("K-S:"),
              uiOutput('ks'),
              h5("A-D:"),
              uiOutput('ad'),
              h5("AIC:"),
              uiOutput('aic'),
              h5("BIC:"),
              uiOutput('bic'),
              h5("Parameters (and their standard errors)"),
              verbatimTextOutput('params')
              
            ),
            mainPanel(
              plotOutput("fitPlot"),
              plotOutput("quantilePlot"),
              plotOutput("finalPlot") 
            )
          )
     ),
  tabPanel("Compare Distributions",
     plotOutput("compare", width = "50%")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    dataInput <- reactive({ 
      data <- read_excel(paste0("data/", input$station, ".xls"), sheet = 1) 
      out <- preprocessing(data)
      data_processed <<- out$processed
      data_sorted <<- out$sorted
      noutliers <<- out$noutliers
      aps <<- out$aps
      return_periods <<- (dim(data_sorted)[1]+1)/data_sorted$rank    # Weibull plotting position   
      if (input$distr == "norm") {
        estim <<- normal(aps, input$method, log = FALSE, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "lognorm") {
        estim <<- normal(log(aps), input$method, log = TRUE, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "expo") {
        estim <<- exponential(aps, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "gam") {
        estim <<- gam(aps, input$method, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "p3") {
        estim <<- p3(aps, input$method, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "lp3") {
        estim <<- lp3(log(aps), input$method, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "gum") {
        estim <<- gumb(aps, input$method, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "wei") {
        estim <<- weibull(aps, input$method, rp = return_periods, plotx = xseq)
      }
      })
    
    output$outliers <- renderText({
      # report number of outliers
      dataInput()
      noutliers
    })
    
    output$timePlot <- renderPlot({
      dataInput()
      if (input$plottype == "Time series") {
          # time series plot
          plot(data_processed$Year, data_processed$aps, type="l", xlab = "Year", ylab = "Maximum discharge (cf/s)", main = "Time series of the raw data")
      }
      
      if (input$plottype == "Density") {  
          # density plot
          hist(data_processed$aps, breaks=15, main = "Distribution of annual peak discharge", xlab = "Annual peak streamflow (cf/s)", prob = TRUE)
          lines(density(data_processed$aps), col = "blue", lwd = 2)
      }
  })

    output$fitPlot <- renderPlot({
      # fitted density plot
      dataInput()
      hist(aps, breaks=15, main = "Distribution of annual peak discharge", xlab = "Annual peak streamflow (cf/s)", prob = TRUE)
      lines(density(aps), col = "blue", lwd = 2)
      lines(xseq, estim$ploty, col = "red" , type="l", lwd = 2)
    })
    
    output$quantilePlot <- renderPlot({
      # q-q plot
      dataInput()
      qqplot(aps, estim$xt, xlab="Observed peaks (cfs)", ylab="Estimated peaks (cfs)", main="Q-Q plot")
      abline(c(0,1), col="red", lty=1, lwd=2)
      legend("topleft", legend=c("1-1 line", "Regression line", " 95% confidence band\n(based on K-S statistic)"),
             col=c("red", "grey", "grey"), lty=c(1,1,2), lwd=c(2,1,1), cex=0.95)
    })
    
    output$finalPlot <- renderPlot({
      # plot of estimated discharge against return period
      dataInput()
      plot(log(return_periods), data_sorted$aps, xlab = "Return period (years)c", ylab = "Discharge (cfs)", 
           main = "Estimated peak discharge vs log return period", xaxt="n")
      lines(log(return_periods), estim$xt, col=2)
      at.x <- log(c(1,5,10,25,50,75,100))
      axis(1, at=at.x, labels=round(exp(at.x), 0), las=1)
    })
    
    output$ks <- renderPrint({
      # Kolmogorov-Smirnov
      dataInput()
      cat(format(signif(suppressWarnings(ks.test(estim$xt, aps)$statistic), 3), nsmall=3))
    })
    
    output$rmse <- renderPrint({
      # root-mean-square-error
      dataInput()
      cat(format(signif(sqrt(sum((sort(aps, decreasing = TRUE)-estim$xt)^2)/length(aps)), 5), nsmall=1))
    })
    
    output$aic <- renderPrint({
      # Akaike Information Criterion
      dataInput()
      if(input$method != "mle"){
        cat("Only available for MLE") 
      } else {
      cat(round(-2*estim$likelihood+2*estim$parno, digits=1))
      }
    })
    
    output$bic <- renderPrint({
      # Bayesian Information Criterion
      dataInput()
      if(input$method != "mle"){
       cat("Only available for MLE")
      } else {
        cat(round(-2*estim$likelihood+estim$parno*log(length(aps)), digits=1))
      }
    })
    
    output$ad <- renderPrint({
      # Anderson-Darling
      dataInput()
      cat(format(signif(estim$ad, 4), nsmall=3))
    })
    
    output$params <- renderPrint({
      # estimated parameters
      dataInput()
      names <- names(estim$par)
      for (i in 1:estim$parno){
        cat(names[i], ": ", format(signif(unlist(estim$par)[i], 5), nsmall=2), " (", format(signif(unlist(estim$sd)[i],4), nsmall=2),")")
        cat("\n")
      }
    })
    
    output$compare <- renderPlot({
      # plot comparing all distributions
      plot_combined(input$station)
    })

}

# Run the application 
shinyApp(ui = ui, server = server)

