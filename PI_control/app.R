# Monte Carlo simulation of Titration Paradox for PI controller
# July, 2020
# St.Gallen, Switzerland
# Thomas Schnider, Prof. Dr. med.
# Charles Minto, MB ChB

# code update March 2024 to demonstrate proportional integral control
#
# Elie Sarraf, M.D. C.M.

list.of.packages <- c("shiny","shinyBS","graphics")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


require(shiny)
require(shinyBS)
require(graphics)
require(zeallot)

# J Pharmacokinet Pharmacodyn. 2010 Jun;37(3):243-55.
# A two-compartment effect site model describes the bispectral index after different rates of propofol infusion.
#
# Björnsson MA1, Norberg A, Kalman S, Karlsson MO, Simonsson US.
# Slope
pGamma <- 2.93
# SD from the Björnsson study (IIV of D50)
mSD <- 0.18
# D50 "normalized"
D50_TV <- 1
#
maxSteps <- 25

Emax <- function(D, D50, g, maxEff = 1, zeroEff = 0) {
  hill <- D^g / (D50^g + D^g)
  return(zeroEff + (maxEff - zeroEff) * hill)
}

# Emax solved for D
revPD <- function(E, E0, Emax, D50, g) {
  return(((E - E0) * D50^g / (Emax - E))^(1 / g))
}
# First Derivative of effect - used for calculation of required change of input
# Thanks to Mathematica: "D[(D/D50)^g/(1 + (D/D50)^g), D]"
d1dD <- pGamma / (4 * D50_TV)

# Log normal variability - for input and between steps
logNormalVar <- function(value, mSD = 0.1, isRV = FALSE) {
  #  rnormA <- repeatable(rnorm)
  p <- value * exp(rnorm(length(value), 0, mSD))
  ifelse(isRV, return(p), return(value))
}


# Calulation of delta input: Proportional change
calc_dose <- function(measEff, currentDose, targetEff, hill1Deriv, startVal) {
  # difference of dose corresp with effect along the linearized Hill curve
  # ... E = (E_50 - E)/d1D1
  error_calc <- (targetEff - measEff) / hill1Deriv
  

  # difference between  current dose and expected dose
  # typical dose for targetEff. targetEff = 0.5 (fix)
  # D_TV <- revPD(targetEff,0,1,D50_TV,pGamma)
  #deviationFromMeanDose_Dose <- D50_TV - currentDose
  
  # Adding integral component instead of modified difference component
  integral_error <- error_calc +startVal$prior_error
  # deltaDoseForNextStep = f(deltaEffect)* scalar1 + f(deltaDose) * scalar2
  #return(list(deviationFromMeanEffect_Dose * startVal$scaleD_Effect + deviationFromMeanDose_Dose * startVal$scaleD_Dose,deviationFromMeanDose_Dose))
  
  return(list(error_calc * startVal$prop_factor + integral_error * startVal$integral_facor,integral_error))
}


# Titrates towards target effect nSteps times
titrate <- function(nSteps, inputSD, targetEffect, D50_ind, startVal) {
  startDoseSD <- startVal$sd
  startDose <- startVal$startDose
  res <- do.call(rbind, apply(
    as.matrix(D50_ind), 1, function(D50) {

      # Empty nSteps by 4 col matrix for result.
      Dose_Effect <- matrix(nrow = nSteps, ncol = 4)
      colnames(Dose_Effect) <- c("Dose", "Effect", "D50", "gamma")
      # Start input
      Dose_Effect[1, 1] <- logNormalVar(startDose, startDoseSD, TRUE)
      #
      D50 <- logNormalVar(D50, mSD = inputSD, isRV = TRUE)
      Dose_Effect[1, 2] <- Emax(Dose_Effect[1, 1], D50, pGamma)
      Dose_Effect[1, 3] <- D50
      Dose_Effect[1, 4] <- pGamma
	  startVal$prior_error = 0
      for (i in 2:nSteps) {
        # variability due to e.g. stimulation
        D50 <- logNormalVar(D50, mSD = inputSD, isRV = TRUE)
        # noise in effect measure or not making use of measurement
        effectToUse <- Dose_Effect[i - 1, "Effect"] + rnorm(length(Dose_Effect[i - 1, "Effect"]), sd = startVal$effectError)

        c(calcDose, temp) %<-% calc_dose(effectToUse, Dose_Effect[i - 1, "Dose"], targetEffect, d1dD, startVal)
        startVal$prior_error <- temp
        newDose <- Dose_Effect[1, "Dose"] + round(calcDose, 20)
        if (newDose < 0) newDose <- 0

        Eff <- Emax(newDose, D50, pGamma)
        Dose_Effect[i, "Dose"] <- newDose
        Dose_Effect[i, "Effect"] <- Eff
        Dose_Effect[i, 3] <- D50
        Dose_Effect[i, 4] <- pGamma
      } # End for
	  

      return(as.data.frame(Dose_Effect[1:nSteps, ]))
    } # End function
  )) # End do.call

  colnames(res) <- c("Dose", "Effect", "D50", "gamma")
  return(res)
}

# Text for ui
textStartDose <- (HTML(paste0("Initial dose as a fraction of D", tags$sub("50"))))

textTitration <- HTML(paste0(
  "For scrolling through the titration steps:<br/>",
  "At each step the new dose is calculated based on the difference between measured ",
  "and target effect and current dose and expected dose:<br/> ",
  "&#916;D = &#946;", tags$sub("1"), "(E", tags$sub("T"), "-E", tags$sub("i"), ")/k",
  "+ &#946;", tags$sub("2"), "(D", tags$sub("T"), "-D", tags$sub("i"), ")<br/>",
  "&#946;", tags$sub("1"), " and &#946;", tags$sub("2"), " are scalars - see below.<br/>",
  "E", tags$sub("T"), " is the target effect <br/> E", tags$sub("i"),
  " is the individual’s current effect<br/>",
  "D", tags$sub("T"), " is the expected or average dose <br/> D", tags$sub("i"),
  " is the individual’s current dose<br/>",
  "k is the first derivative (slope) of the Emax function evaluated at D",
  tags$sub("50"), " (red line)"
))

textEffect <- HTML(paste0(
  "Scalar &#946;",
  tags$sub("1"), " in:                    <br/>",
  "&#916;D = &#946;", tags$sub("1"), "(E", tags$sub("T"), "-E", tags$sub("i"), ")/k",
  "+ &#946;", tags$sub("2"), "(D", tags$sub("T"), "-D", tags$sub("i"), ")"
))
textDose <- paste0(
  "Scalar &#946;",
  tags$sub("2"), " in:                    <br/>",
  "&#916;D = &#946;", tags$sub("1"), "(E", tags$sub("T"), "-E", tags$sub("i"), ")/k",
  "+ &#946;", tags$sub("2"), "(D", tags$sub("T"), "-D", tags$sub("i"), ")<br/>",
  "This reflects our tendency in titration to e.g. reduce the dose more ",
  "agressively if effect is high AND the current dose higher than expected"
)
textEffectError <- paste0(
  "Variability or error in the measurement. This is technically the same as not relying ",
  "on the measured effect for dosing. The bigger SD is, the more random is the dose i.e. ",
  "independent of the measurement"
)

# ui:

shinyApp(
  ui <- fluidPage(

    # Application title
    titlePanel("TitrationParadoxR"), 

    # Sidebar with a slider input for number of bins
    sidebarLayout(
      sidebarPanel(
        sliderInput("stepToPrint",
          "Step to print:",
          min = 0,
          max = maxSteps,
          value = 0
        ),
        bsTooltip("stepToPrint", textTitration,
          options = list(container = "body")
        ),
        sliderInput("doseScalarEffect",
          HTML(paste0(
            "Scalar ", "&#946;",
            tags$sub("1"), " (scales ", "P)"
          )),
          min = 0,
		  # max changed from 1 to 2
          max = 2,
          step = .01,
          value = 0.6
        ),
        bsTooltip("doseScalarEffect", textEffect,
          options = list(container = "body")
        ),
        sliderInput("doseScalarDose",
          HTML(paste0(
            "Scalar ", "&#946;",
            tags$sub("2"), " (scales ", "I)"
          )),
          min = 0,
          max = 2, # change range to accomadate integral effect
          step = 0.01,
          value = 0.6
        ),
        bsTooltip("doseScalarDose", textDose,
          options = list(container = "body")
        ),

        sliderInput("inputSD",
          "Between steps variability (SD):",
          min = 0.0,
          max = 0.1,
          step = 0.001,
          value = .015
        ),
        sliderInput("effectError",
          "Inaccuracy or ignorance of measurement (SD):",
          min = 0.0,
          max = 0.3,
          step = 0.01,
          value = 0
        ),
        bsTooltip("effectError", textEffectError,
          options = list(container = "body")
        ),

        sliderInput("startDose",
          textStartDose,
          min = 0.05,
          max = 2,
          step = .05,
          value = 1
        ),
        bsTooltip("startDose", textStartDose, options = list(container = "body")),
        sliderInput("startSDDose",
          "Variability of initial dose (SD):",
          min = 0,
          max = 0.4,
          step = 0.01,
          value = 0
        ),

        sliderInput("numSim",
          "Number of patients:",
          min = 100,
          max = 1000,
          step = 100,
          value = 500
        )
      ),


      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("effDosePlot", width = "100%", height = "500px")
      )
    )
  ),

  # server:

  server <- function(input, output) {
    D50_ind <- reactive(
      return(D50_TV * exp(rnorm(as.numeric(input$numSim), 0, mSD)))
    )


    input_Changed <- reactive({
      # tempScaler <- 0
      # if (input$acknowledgeDose) tempScaler <- 0.2
      #


      startVal <- list(
        sd = input$startSDDose,
        # startDose = input$startDose, scaleD_Dose = tempScaler,
        startDose = input$startDose, integral_facor = input$doseScalarDose,
        prop_factor = input$doseScalarEffect, effectError = input$effectError, prior_error = 0
      )
      nSubj <- as.numeric(input$numSim)

      nSteps <- maxSteps + 1

      inputSD <- as.numeric(input$inputSD)
      targetEff <- 0.5
      D50_ind <- D50_TV * exp(rnorm(nSubj, 0, mSD))
      list(
        retVals = list(startVal = startVal, nSteps = nSteps, nSubj = nSubj, D50_ind = D50_ind(), targetEff = targetEff, inputSD = inputSD),
        Dose_EffectVals = titrate(nSteps, inputSD, targetEff, D50_ind(), startVal)
      )
    })

    output$effDosePlot <- renderPlot({
      if (is.null(input$stepToPrint)) {
        return()
      }
      idx <- as.numeric(input$stepToPrint) + 1


      Dose_EffectValues <- input_Changed()$Dose_EffectVals
      nSteps <- input_Changed()$retVals$nSteps
      nSubj <- input_Changed()$retVals$nSubj
      targetEff <- input_Changed()$retVals$targetEff
      D50_ind <- input_Changed()$retVals$D50_ind
      inputSD <- input_Changed()$retVals$inputSD

      par(mfrow  = c(2,1),mai = c(1, 1, 0, 0.42), cex = 1.5)
      #par()
      plot(0, 0, xlab = "", ylab = "", xlim = c(0, 2), ylim = c(0, 1), axes = F, type = "n")
      axis(1, at = seq(0, 2, 0.1))
      axis(2, at = seq(0, 1, 0.1))
      mtext(bquote("Dose (fraction of D"[50] ~ ")"), 1, 2, cex = 1.5)
      mtext("Effect", 2, 2, cex = 1.5)

      # is the selected step!
      stepEffect <- Dose_EffectValues$Effect[seq(idx, nSteps * nSubj, nSteps)]
      stepDose <- Dose_EffectValues$Dose[seq(idx, nSteps * nSubj, nSteps)]

      points(stepDose, stepEffect)
      # linear regression for this step

      # For internal use:
      # writeSimuData(idx,input,stepEffect,stepDose)

      if (idx > 1) {
        mLm <- lm(stepEffect ~ stepDose)
        abline(mLm, lwd = 1, col = "blue")
      }

      
      # Plot the PD curve
      Effs <- seq(0.01, 0.99, 0.01)
      Doss <- seq(0.01, 1.99, 0.01)

      # Solved Hill Emax
      lines(Doss, Emax(Doss, D50_TV, pGamma), lwd = 3, lty = 3)
      lines(c(1, 1), c(0, 1), col = "red", lwd = 4, lty = 4)
      # Horizontal line of target effect
      lines(c(0, 2), c(input_Changed()$retVals$targetEff, input_Changed()$retVals$targetEff), col = "red", lwd = 4, lty = 4)


      qEff <- round(quantile(stepEffect, c(0.005, 0.995)), 2)
      text(.4, .75, paste("99% effect between\n ", qEff[1], "and", qEff[2]), cex = .75)
      text(1.8, .3, paste("CVe:\n ", round(100 * sd(stepEffect) / mean(stepEffect), 2), "%"), cex = .75)
      text(1.8, .15, paste("CVd:\n ", round(100 * sd(stepDose) / mean(stepDose), 2), "%"), cex = .75)
      # Line with first derivative
      # (.5 - d1dD) + Dose * d1dD
      y1 <- (.5 - d1dD) + 0.3 * d1dD
      y2 <- (.5 - d1dD) + 1.8 * d1dD
      lines(c(0.3, 1.8), c(y1, y2), lwd = 2, col = "red")
      
      # ES addition to code
      # demonstrate and plot the error function and add to plot
      corr_mat <- vector(length = nSteps)
      for (idx1 in 2:nSteps){
        stepEffect1 <- Dose_EffectValues$Effect[seq(idx1, nSteps * nSubj, nSteps)]
        stepDose1 <- Dose_EffectValues$Dose[seq(idx1, nSteps * nSubj, nSteps)]
        corr_mat[idx1]= cor(x=stepDose1,y=stepEffect1)
      }
      print(sum(corr_mat))
      plot(1:(nSteps-1), corr_mat[2:nSteps], xlab = "", ylab = "", xlim = c(0, nSteps), ylim = c(-1, 1), axes = F, type = "p")
      axis(1, at = seq(0, 25, 1))
      axis(2, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
      mtext("Correlation", 2, 2, cex = 1.5)
      mtext(bquote("Steps"), 1, 2, cex = 1.5)
      abline(h = 0, col="red", lwd=3, lty=2)
      #print(corr_mat)

    })
  }
)
