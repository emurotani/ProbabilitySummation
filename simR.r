#
# An implementation of the OPI that simulates responses using 
# Cummulative Gaussian variability with supplied standard deviation.
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: June 2012
#
# Modified 20 Jul 2014: added maxStim argument for cdTodB conversion
#
# Copyright 2012 Andrew Turpin
# This program is part of the OPI (http://perimetry.org/OPI).
# OPI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#' @rdname opiClose
#' @return
#' \subsection{simRaussian}{
#'   DETAILS
#' }
#' @examples
#' chooseOpi("simRaussian")
#' if (!is.null(opiInitialize(sd=2)))
#'   stop("opiInitialize failed")
#' if (!is.null(opiClose()))
#'   stop("opiClose failed, which is very surprising!")
simR.opiClose         <- function() { return(NULL) }

#' @rdname opiQueryDevice
#' @examples
#' chooseOpi("simRaussian")
#' if (!is.null(opiInitialize(sd=2)))
#'   stop("opiInitialize failed")
#' print(opiQueryDevice())
simR.opiQueryDevice   <- function() { return (list(type="simRaussian", isSim=TRUE)) }

if (exists(".OpiEnv") && !exists("simR", where=.OpiEnv))
    assign("simR", new.env(2), envir=.OpiEnv)

################################################################################
# Input
#   sd standard deviation for the Gaussian
#   display Dimensions of plot area to display stim. No display if NULL
#
# Return NULL if succesful, string error message otherwise  
################################################################################
#' @rdname opiInitialize
#' @param sd standard deviation for the Gaussian
#' @details
#' \subsection{simRaussian}{
#'   \code{opiInitialize(sd, display=NA, maxStim=10000/pi)}
#'   
#'   If the chosen OPI implementation is \code{simRaussian}, then \code{sd} is the
#'   standard deviation value that the simulator will use for the slope/spread of
#'   the psychometric function.
#'   
#'   \code{display} and \code{maxStim} is as for SimHenson.
#' }
#' @examples
#' # Set up a simulation using a psychometric function that is
#' # a cumulative gaussian of standard deviation 2
#' chooseOpi("simRaussian")
#' if (!is.null(opiInitialize(sd=2)))
#'   stop("opiInitialize failed")
simR.opiInitialize <- function(sd = 2, display = NA, maxStim = 10000 / pi) {
    if (!is.numeric(sd) || (sd < 0)) {
        msg <- paste("Invalid standard deviation in opiInitialize for simRaussian:",sd)
        warning(msg)
        return(msg)
    }

    df.est.par <- read.csv("3532912_estimated_parameters.csv")

    mu <- 18 # <0dB
    mu <- c(mu, df.est.par[1:36, 3]) # 0~35dB
    .OpiEnv$simR$mu <- mu

    sigma <- 0.50
    sigma <- c(sigma, df.est.par[37:72, 3])
    .OpiEnv$simR$sigma <- sigma
    
    lamda <- 1.0
    lamda <- c(lamda, df.est.par[73:108, 3])
    .OpiEnv$simR$lamda <- lamda
  
    .OpiEnv$simR$sd <- sd
    .OpiEnv$simR$maxStim <- maxStim

    if (simDisplay.setupDisplay(display))
        warning("opiInitialize (simRaussian): perhaps display parameter does not contain 4 numbers?")

    return(NULL)
}

################################################################################
# Set background of plot area to col, color of gird lines to gridCol
################################################################################
#' @rdname opiSetBackground
#' @details
#' \subsection{simRaussian}{
#'   \code{opiSetBackground(col, gridCol)}
#'   
#'   \code{col} is the background color of the plot area used for displaying
#'   stimuli, and \code{gridCol} the color of the gridlines. Note the plot area
#'   will only be displayed if \code{opiInitialize} is called with a valid display
#'   argument.
#' }
#' @examples
#' chooseOpi("simRaussian")
#' if (!is.null(opiInitialize(sd=2)))
#'   stop("opiInitialize failed")
#' if (!is.null(opiSetBackground(col="white",gridCol="grey")))
#'   stop("opiSetBackground failed, which is very surprising!")
simR.opiSetBackground <- function(col, gridCol) { 
    simDisplay.setBackground(col, gridCol)
    return(NULL) 
}

#' @rdname opiPresent
#' @details
#' \subsection{simRaussian}{
#'   \code{opiPresent(stim, nextStim=NULL, fpr=0.03, fnr=0.01, tt=30)}
#'   
#'   If the chosen OPI implementation is \code{simRaussian}, then the response
#'   to a stimuli is determined by sampling from a Frequency-of-Seeing (FoS)
#'   curve (also known as the psychometric function) with formula
#'   \code{fpr+(1-fpr-fnr)*(1-pnorm(x, tt, simR.global.sd))}, where \code{x}
#'   is the stimulus value in Humphrey dB, and \code{simR.global.sd} is
#'   set with \code{opiInitialize}.
#' }
simR.opiPresent <- function(stim, nextStim=NULL, fpr=0.03, fnr=0.01, tt=30) { UseMethod("simR.opiPresent") }
setGeneric("simR.opiPresent")

simR.opiPresent.opiStaticStimulus <- function(stim, nextStim=NULL, fpr=0.03, fnr=0.01, tt=30) {
    if (!exists("sd", envir=.OpiEnv$simR)) {
        return ( list(
            err = "opiInitialize(sd) was not called before opiPresent()",
            seen= NA,
            time= NA 
        ))
    }

    if (is.null(stim))
        stop("stim is NULL in call to opiPresent (using simRaussian, opiStaticStimulus)")

    if (as.numeric(.OpiEnv$simR$sd) <= 0)
        warning("sd is <= 0 in simRaussian call to opiPresent")
    
    # print(.OpiEnv$simR$mu[tt])
    
    A <- fpr
    B <- 1-fpr-fnr-as.numeric(.OpiEnv$simR$lamda[tt+2])
    C <- cdTodb(stim$level, .OpiEnv$simR$maxStim)
    D <- as.numeric(.OpiEnv$simR$mu[tt+2])
    E <- as.numeric(.OpiEnv$simR$sigma[tt+2])

    prSeeing <- A + B*(1-pnorm(C, mean=D, sd=E))

    # prSeeing <- fpr + (1-fpr-fnr)*(1-pnorm(cdTodb(stim$level, .OpiEnv$simR$maxStim), mean=tt, sd=as.numeric(.OpiEnv$simR$sd)))

    simDisplay.present(stim$x, stim$y, stim$color, stim$duration, stim$responseWindow)

    return ( list(
        err = NULL,
        seen= runif(1) < prSeeing,
        time= 0
    ))
}#

########################################## TO DO
simR.opiPresent.opiTemporalStimulus <- function(stim, nextStim=NULL, ...) {
    stop("ERROR: haven't written simR temporal persenter yet")
}

simR.opiPresent.opiKineticStimulus <- function(stim, nextStim=NULL, ...) {
    stop("ERROR: haven't written simR kinetic persenter yet")
}
