isTRUEorFALSE = function(x) {
  isTRUE(x) || isFALSE(x)
}

as.susy = function(x) {
  if (inherits(x, "susy"))
    return(x)
  if (!is.list(x))
    stop("only list class objects can be turned into susy objects")
  class(x) = unique(c("susy", class(x)))
  x
}

## loops using internal crosscor function over column pairs
susy = function(x, segment, Hz, maxlag=3L, permutation=FALSE, restrict.surrogates=FALSE, surrogates.total=500) {
  if (!is.data.frame(x))
    stop("'x' must be a data.frame")
  if (is.null(names(x)))
    stop("'x' must have named columns")
  if (!all(vapply(x, is.numeric, FALSE, USE.NAMES=FALSE)))
    stop("'x' must have numeric columns")
  if (!isTRUEorFALSE(permutation))
    stop("'permutation' must be TRUE or FALSE")
  if (!is.numeric(segment) || length(segment)!=1L || is.na(segment))
    stop("'segment' must be scalar non-NA numeric")
  if (!is.numeric(Hz) || length(Hz)!=1L || is.na(Hz))
    stop("'Hz' must be scalar non-NA numeric")
  if (!is.numeric(maxlag) || length(maxlag)!=1L || is.na(maxlag))
    stop("'maxlag' must be scalar non-NA numeric")

  nx = ncol(x)
  if (nx < 2L)
    stop("'x' must have at least 2 columns")
  if (!permutation && nx%%2L)
    stop("When 'permutation' is FALSE then 'x' must have even number of columns")

  segment = as.integer(segment)
  Hz = as.integer(Hz)
  maxlag = as.integer(maxlag)

  if (segment < 2L*maxlag)
    stop("'segment' must not be smaller than '2 * maxlag'")
  if (segment > nrow(x)/2)
    stop("'segment' must not be greater than 'nrow(x)/2'")

  maxlagHz = maxlag*Hz
  lagtimes2 = maxlagHz*2L + 1L
  range = segment * Hz

  if (!permutation) {
    pairs = matrix(seq_len(nx), ncol=2L, byrow=TRUE)
  } else {
    pairs = t(combn(seq_len(nx), 2L))
  }

  crosscor = function(ipair) {
    vars = names(x)[pairs[ipair,]]
    a = x[, pairs[ipair,1L]]
    b = x[, pairs[ipair,2L]]
    cc = complete.cases(a, b)
    a = a[cc]
    b = b[cc]
    if (!length(a))
      stop("empty data after filtering out missing values")
    size = length(a) ## same as length of b
    numberEpochen = round(size/range-0.499999)
    if (restrict.surrogates) {
      anzahlPseudosProEpoche = floor(surrogates.total/numberEpochen)
      if (anzahlPseudosProEpoche < 1L) {
        warning("Number of pseudos per total is bigger than number number of segments (length/range-0.5), setting pseudos per segment to 1")
        anzahlPseudosProEpoche = 1L
      }
      if (anzahlPseudosProEpoche >= numberEpochen) {
        warning("surrogates.total is bigger than number of segments (length/range-0.5), setting pseudos per segment to number of segments - 1")
        anzahlPseudosProEpoche = numberEpochen-1L
      }
      if (anzahlPseudosProEpoche == numberEpochen-1L) {
        restrict.surrogates = FALSE
      }
    } else {
      anzahlPseudosProEpoche = numberEpochen-1L
    }

    eps = array(0, dim=c(2L, numberEpochen, range))
    i2 = 1L
    for (i0 in 1:numberEpochen) {
      for (i1 in 1:range) {
        eps[1L,i0,i1] = a[i2]
        eps[2L,i0,i1] = b[i2]
        i2 = i2 + 1L
      }
    }

    meanccorrReal = meanccorrPseudo = meanccorrRealZ = meanccorrRealZNotAbs = meanccorrPseudoZ = meanccorrPseudoZNotAbs =
      vector("double", lagtimes2)
    nReal = nPseudo =
      vector("double", numberEpochen)
    if (restrict.surrogates) {
      for (i in seq_len(numberEpochen)) {
        result = ccf(eps[1L,i,], eps[2L,i,], lag.max=maxlagHz, plot=FALSE)
        if (!is.nan(result$acf[1L])) {
          for (j in seq_len(lagtimes2)) {
            meanccorrRealZ[j] = meanccorrRealZ[j] + abs(0.5*log((1+result$acf[j])/(1-result$acf[j])))
            meanccorrRealZNotAbs[j] = meanccorrRealZNotAbs[j] + 0.5*log((1+result$acf[j])/(1-result$acf[j]))
            meanccorrReal[j] = meanccorrReal[j] + abs(result$acf[j])
            nReal[i] = nReal[i] + abs(result$acf[j])
          }
          nReal[i] = nReal[i]/lagtimes2
        }
      }
      isOcc = vector("integer", numberEpochen)
      n = 0L
      v = 1L
      for (i in seq_len(numberEpochen)) {
        for (h in seq_len(anzahlPseudosProEpoche)) {
          repeat {
            v = floor(runif(1L, 1L, numberEpochen+1L))
            if (isOcc[v] != i && v != i){
              isOcc[v] = i
              break
            }
          }
          n = n + 1L
          result = ccf(eps[1L,i,], eps[2L,v,], lag.max=maxlagHz, plot=FALSE)
          if (!is.nan(result$acf[1L])) {
            for (j in seq_len(lagtimes2)) {
              meanccorrPseudoZ[j] = meanccorrPseudoZ[j] + abs(0.5*log((1+result$acf[j])/(1-result$acf[j])))
              meanccorrPseudoZNotAbs[j] = meanccorrPseudoZNotAbs[j] + 0.5*log((1+result$acf[j])/(1-result$acf[j]))
              meanccorrPseudo[j] = meanccorrPseudo[j] + abs(result$acf[j])
              nPseudo[i] = nPseudo[i] + abs(result$acf[j])
            }
          }
        }
        nPseudo[i] = nPseudo[i]/(anzahlPseudosProEpoche*lagtimes2)
      }
      meanccorrPseudo = meanccorrPseudo/n
      meanccorrPseudoZ = meanccorrPseudoZ/n
      meanccorrPseudoZNotAbs = meanccorrPseudoZNotAbs/n
    } else {
      for (i in seq_len(numberEpochen)) {
        for (h in seq_len(numberEpochen)) {
          result = ccf(eps[1L,i,], eps[2L,h,], lag.max=maxlagHz, plot=FALSE)
          if (is.nan(result$acf[1L]))
            next
          for (j in seq_len(lagtimes2)) {
            if (i == h) {
              meanccorrRealZ[j] = meanccorrRealZ[j] + abs(0.5*log((1+result$acf[j])/(1-result$acf[j])))
              meanccorrRealZNotAbs[j] = meanccorrRealZNotAbs[j] + 0.5*log((1+result$acf[j])/(1-result$acf[j]))
              meanccorrReal[j] = meanccorrReal[j] + abs(result$acf[j])
              nReal[i] = nReal[i] + abs(result$acf[j])
            } else {
              meanccorrPseudoZ[j] = meanccorrPseudoZ[j] + abs(0.5*log((1+result$acf[j])/(1-result$acf[j])))
              meanccorrPseudoZNotAbs[j] = meanccorrPseudoZNotAbs[j] + 0.5*log((1+result$acf[j])/(1-result$acf[j]))
              meanccorrPseudo[j] = meanccorrPseudo[j] + abs(result$acf[j])
              nPseudo[i] = nPseudo[i] + abs(result$acf[j])
            }
          }
        }
        nReal[i] = nReal[i]/lagtimes2
        nPseudo[i] = nPseudo[i]/((numberEpochen-1L)*lagtimes2)
      }
      meanccorrPseudo = meanccorrPseudo/((numberEpochen-1L)*numberEpochen)
      meanccorrPseudoZ = meanccorrPseudoZ/((numberEpochen-1L)*numberEpochen)
      meanccorrPseudoZNotAbs = meanccorrPseudoZNotAbs/((numberEpochen-1L)*numberEpochen)
    }
    meanccorrReal = meanccorrReal/numberEpochen
    meanccorrRealZ = meanccorrRealZ/numberEpochen
    meanccorrRealZNotAbs = meanccorrRealZNotAbs/numberEpochen

    k1 = k1NotAbs = 0L
    for (t in seq_len(lagtimes2)) {
      if (meanccorrRealZ[t] > meanccorrPseudoZ[t])
        k1 = k1 + 1L
      if (meanccorrRealZNotAbs[t] > meanccorrPseudoZNotAbs[t])
        k1NotAbs = k1NotAbs + 1L
    }

    ## results
    list(
      data = setNames(data.frame(a, b), vars),
      lagtimes2.data = data.frame(
        meanccorrReal=meanccorrReal,
        meanccorrPseudo=meanccorrPseudo,
        meanccorrRealZ=meanccorrRealZ,
        meanccorrRealZNotAbs=meanccorrRealZNotAbs,
        meanccorrPseudoZ=meanccorrPseudoZ,
        meanccorrPseudoZNotAbs=meanccorrPseudoZNotAbs
      ),
      segment.data = data.frame(
        nReal=nReal,
        nPseudo=nPseudo
      ),
      params = list(
        segment=segment, numberEpochen=numberEpochen, size=size, maxlag=maxlag, maxlagHz=maxlagHz, lagtimes2=lagtimes2,
        anzahlPseudosProEpoche=anzahlPseudosProEpoche,
        k1=k1, k1NotAbs=k1NotAbs
      )
    )
  }
  ipairs = seq_len(nrow(pairs))
  pair.names = function(i) paste(names(x)[pairs[i,]], collapse="-")
  names(ipairs) = vapply(ipairs, pair.names, "", USE.NAMES=FALSE)
  ans = lapply(ipairs, crosscor)
  as.susy(ans)
}

## this re-uses list subset so we can nicely subset susy object
"[.susy" = function(x, ...) {
  as.susy(NextMethod(x, ...))
}

## plot various types of plots using susy object
plot.susy = function(x, type=c(4, 5), ...) {
  if (!inherits(x, "susy"))
    stop("'x' must be an object of class 'susy'")
  if (!is.numeric(type))
    stop("'type' argument must be numeric")
  type = as.integer(type)
  if (anyNA(type) || any(type > 5L) || any(type < 1L))
    stop("'type' values must be in range of 1 to 5")
  if (anyDuplicated(type))
    stop("'type' values must be unique")

  first.plot = TRUE
  plot1 = function(x) {
    ## take out from susy object so no need to prefix in all code below
    a = x$data[[1L]]
    b = x$data[[2L]]
    variablenname1 = names(x$data)[1L]
    variablenname2 = names(x$data)[2L]
    size = x$params$size
    segment = x$params$segment
    anzahlPseudosProEpoche = x$params$anzahlPseudosProEpoche
    numberEpochen = x$params$numberEpochen
    maxlagHz = x$params$maxlagHz
    meanccorrReal = x$lagtimes2.data$meanccorrReal
    meanccorrPseudo = x$lagtimes2.data$meanccorrPseudo
    meanccorrRealZ = x$lagtimes2.data$meanccorrRealZ
    meanccorrPseudoZ = x$lagtimes2.data$meanccorrPseudoZ
    meanccorrRealZNotAbs = x$lagtimes2.data$meanccorrRealZNotAbs
    meanccorrPseudoZNotAbs = x$lagtimes2.data$meanccorrPseudoZNotAbs
    nReal = x$segment.data$nReal
    nPseudo = x$segment.data$nPseudo
    for (t in type) {
      if (!first.plot) dev.new() else first.plot <<- FALSE
      if (t == 1L) { ## synchronie meanccorrReal meanccorrPseudo
        min0 = min(meanccorrReal, meanccorrPseudo)
        max0 = max(meanccorrReal, meanccorrPseudo)
        title0 = sprintf("Synchrony %s-%s segment: %ss; %s pseudos", variablenname1, variablenname2, segment, anzahlPseudosProEpoche*numberEpochen)
        plot(seq(from = -maxlagHz , to = maxlagHz), meanccorrReal, ylim=c(min0, max0+ 0.19*(max0-min0)), main=title0,
             xlab="lag", ylab="correlation", type="l", col="green", lwd = 4)
        lines(seq(from = -maxlagHz , to = maxlagHz), meanccorrPseudo, col="red", lwd = 4)
        legend("topright", pch = c(3), col = c("green", "red"),legend = c("meanccorr", "meanccorr pseudo"))
      } else if (t == 2L) { ## Epochensynchronie nReal nPseudo
        title0 = sprintf("Segment synchronies %s-%s segment: %ss; %s pseudos", variablenname1, variablenname2, segment, anzahlPseudosProEpoche*numberEpochen)
        plot(seq(from = 1 , to = numberEpochen), nReal, ylim=c(min(nReal, nPseudo), max(nReal, nPseudo)), main=title0,
             xlab="segment", ylab="correlation", type="l", col="green")
        lines(seq(from = 1, to = numberEpochen ), nPseudo, col="red", lwd = 2)
        legend("topright", pch = c(3), col = c("green", "red"),legend = c("cor real", "cor pseudo"))
      } else if (t == 3L) { ## Synchrony meanccorrRealZ meanccorrPseudoZ
        min0 = min(meanccorrRealZ, meanccorrPseudoZ)
        max0 = max(meanccorrRealZ, meanccorrPseudoZ)
        title0 = sprintf("Z-Synchrony %s-%s segment: %ss; %s pseudos", variablenname1, variablenname2, segment, anzahlPseudosProEpoche*numberEpochen)
        plot(seq(from = -maxlagHz , to = maxlagHz), meanccorrRealZ, ylim=c(min0, max0+ 0.19*(max0-min0)),
             main=title0, xlab="lag", ylab="correlation", type="l", col="green", lwd = 4)
        lines(seq(from = -maxlagHz , to = maxlagHz), meanccorrPseudoZ, col="red", lwd = 4)
        legend("topright", pch = c(3), col = c("green", "red"),legend = c("Z(meanccorr)", "Z(meanccorr pseudo)"))
      } else if (t == 4L) { ## Zeige die Verteilung der data an
        title0 = sprintf("Time series %s-%s segment: %ss; %s pseudos", variablenname1, variablenname2, segment, anzahlPseudosProEpoche*numberEpochen)
        plot(seq(from = 1 , to = size), a,ylim=c(min(a, b), max(a, b)), main=title0,
             xlab="Zeit", ylab="Wert", type="l", pch=20, col="green")
        points(seq(from = 1 , to = size), b, type="l", pch=20, col="red", lwd = 2)
        legend("topright", pch = c(20,20), col = c("green", "red"),legend=c(variablenname1, variablenname2))
      } else if (t == 5L) { ## Synchrony meanccorrRealZNotAbs meanccorrPseudoZNotAbs
        min0 = min(meanccorrRealZNotAbs, meanccorrPseudoZNotAbs)
        max0 = max(meanccorrRealZNotAbs, meanccorrPseudoZNotAbs)
        title0 = sprintf("Z-Synchrony not ABS %s-%s segment: %ss; %s pseudos", variablenname1, variablenname2, segment, anzahlPseudosProEpoche*numberEpochen)
        plot(seq(from = -maxlagHz , to = maxlagHz), meanccorrRealZNotAbs, ylim=c(min0, max0+ 0.19*(max0-min0)),
             main=title0, xlab="lag", ylab="correlation", type="l", col="green", lwd = 4)
        lines(seq(from = -maxlagHz , to = maxlagHz), meanccorrPseudoZNotAbs, col="red", lwd = 4)
        legend("topright", pch = c(3), col = c("green", "red"),legend = c("Z(meanccorr not ABS)", "Z(meanccorr pseudo not ABS)"))
      }
    }
  }
  lapply(x, plot1)
  invisible()
}

as.data.frame.susy = function(x, row.names=NULL, optional=FALSE, corr.no.abs=TRUE, ...) {
  if (!inherits(x, "susy"))
    stop("'x' must be an object of class 'susy'")
  if (!isTRUEorFALSE(corr.no.abs))
    stop("'corr.no.abs' must be TRUE or FALSE")

  cols = c("Var1","Var2","n(data)","Z","Z-Pseudo","SD(Z)","SD(Z-Pseudo)","n(lags)","%>Pseudo","n(Segmente)","ES","Z(lead1)","Z(lead2)","ES(lead1)","ES(lead2)")
  if (corr.no.abs)
    cols = c(cols, c("meanZ(in-phase)","meanZ(anti-phase)","Anzahl(in-phase)","Anzahl(anti-phase)","Z(noAbs)","Z(Pseudo-noAbs)","%>Pseudo(noAbs)","ES(noAbs)"))
  df1 = function(x) {
    a = x$data[[1L]]
    b = x$data[[2L]]
    variablenname1 = names(x$data)[1L]
    variablenname2 = names(x$data)[2L]
    size = x$params$size
    segment = x$params$segment
    anzahlPseudosProEpoche = x$params$anzahlPseudosProEpoche
    numberEpochen = x$params$numberEpochen
    maxlagHz = x$params$maxlagHz
    k1 = x$params$k1
    k1NotAbs = x$params$k1NotAbs
    lagtimes2 = x$params$lagtimes2
    meanccorrReal = x$lagtimes2.data$meanccorrReal
    meanccorrPseudo = x$lagtimes2.data$meanccorrPseudo
    meanccorrRealZ = x$lagtimes2.data$meanccorrRealZ
    meanccorrPseudoZ = x$lagtimes2.data$meanccorrPseudoZ
    meanccorrRealZNotAbs = x$lagtimes2.data$meanccorrRealZNotAbs
    meanccorrPseudoZNotAbs = x$lagtimes2.data$meanccorrPseudoZNotAbs
    nReal = x$segment.data$nReal
    nPseudo = x$segment.data$nPseudos
    l = list(
      variablenname1, variablenname2, size, mean(meanccorrRealZ), mean(meanccorrPseudoZ), sd(meanccorrRealZ), sd(meanccorrPseudoZ), length(meanccorrRealZ),
      100*k1/(length(meanccorrRealZ)),
      numberEpochen,
      (mean(meanccorrRealZ)-mean(meanccorrPseudoZ))/sd(meanccorrPseudoZ),
      mean(meanccorrRealZ[1:maxlagHz]),
      mean(meanccorrRealZ[(maxlagHz+2):lagtimes2]),
      (mean(meanccorrRealZ[1:maxlagHz])-mean(meanccorrPseudoZ[1:maxlagHz]))/sd(meanccorrPseudoZ[1:maxlagHz]),
      (mean(meanccorrRealZ[(maxlagHz+2):lagtimes2])-mean(meanccorrPseudoZ[(maxlagHz+2):lagtimes2]))/sd(meanccorrPseudoZ[(maxlagHz+2):lagtimes2])
    )
    if (corr.no.abs) {
      l = c(l, list(
        mean(meanccorrRealZNotAbs[which(meanccorrRealZNotAbs >= 0)]),
        mean(meanccorrRealZNotAbs[which(meanccorrRealZNotAbs < 0)]),
        length(meanccorrRealZNotAbs[which(meanccorrRealZNotAbs >= 0)]),
        length(meanccorrRealZNotAbs[which(meanccorrRealZNotAbs < 0)]),
        mean(meanccorrRealZNotAbs),
        mean(meanccorrPseudoZNotAbs),
        100*k1NotAbs/(length(meanccorrRealZNotAbs)),
        (mean(meanccorrRealZNotAbs)-mean(meanccorrPseudoZNotAbs))/sd(meanccorrPseudoZNotAbs)
      ))
    }
    setNames(as.data.frame(l), cols)
  }
  do.call("rbind", lapply(unname(x), df1))
}

print.susy = function(x, corr.no.abs=TRUE, legacy=FALSE, ...) {
  if (!inherits(x, "susy"))
    stop("'x' must be an object of class 'susy'")
  if (!isTRUEorFALSE(corr.no.abs))
    stop("'corr.no.abs' must be TRUE or FALSE")
  if (!isTRUEorFALSE(legacy))
    stop("'legacy' must be TRUE or FALSE")

  if (!legacy) {
    df = as.data.frame(x, corr.no.abs=corr.no.abs)
    print(df, ...)
  } else {
    first.print = TRUE
    print1 = function(x) {
      ## take out from susy object so no need to prefix in all code below
      a = x$data[[1L]]
      b = x$data[[2L]]
      variablenname1 = names(x$data)[1L]
      variablenname2 = names(x$data)[2L]
      size = x$params$size
      segment = x$params$segment
      anzahlPseudosProEpoche = x$params$anzahlPseudosProEpoche
      numberEpochen = x$params$numberEpochen
      maxlagHz = x$params$maxlagHz
      k1 = x$params$k1
      k1NotAbs = x$params$k1NotAbs
      lagtimes2 = x$params$lagtimes2
      meanccorrReal = x$lagtimes2.data$meanccorrReal
      meanccorrPseudo = x$lagtimes2.data$meanccorrPseudo
      meanccorrRealZ = x$lagtimes2.data$meanccorrRealZ
      meanccorrPseudoZ = x$lagtimes2.data$meanccorrPseudoZ
      meanccorrRealZNotAbs = x$lagtimes2.data$meanccorrRealZNotAbs
      meanccorrPseudoZNotAbs = x$lagtimes2.data$meanccorrPseudoZNotAbs
      nReal = x$segment.data$nReal
      nPseudo = x$segment.data$nPseudos
      if (first.print) {
        cat("Var1 ")
        cat("Var2 ")
        cat("n(data) ")
        cat("Z ")
        cat("Z-Pseudo ")
        cat("SD(Z) ")
        cat("SD(Z-Pseudo) ")
        cat("n(lags) ")
        cat("%>Pseudo ")
        cat("n(Segmente) ") #10
        # cat("n(Pseudo) ")
        cat("ES ")
        cat("Z(lead1) ")
        cat("Z(lead2) ") #14
        # cat("mean(Z-Pseudo)<lag0 ")
        # cat("mean(Z-Pseudo)>lag0 ") #16
        cat("ES(lead1) ") #CCF wenn der erste (spalte 1) Zahlenstrang fuhrend ist
        cat("ES(lead2) ") #CCF wenn der zweite (spalte 2) Zahlenstrang fuhrend ist
        if (corr.no.abs) {
          cat("meanZ(in-phase) ") #durchschnitt aller positiven korr
          cat("meanZ(anti-phase) ") #durchschnitt aller negativen korr
          cat("Anzahl(in-phase) ") #anzahl aller positiven korr
          cat("Anzahl(anti-phase) ") #anzahl aller negativen korr
          cat("Z(noAbs) ") #25
          cat("Z(Pseudo-noAbs) ")
          cat("%>Pseudo(noAbs) ") # % meanccorrRealZNotAbs die grosser als meanccorrRealZNotAbs sind
          cat("ES(noAbs) ") #effektstarke
        }
        cat("\n")
        first.print <<- FALSE
      }
      cat(variablenname1,"")
      cat(variablenname2,"")
      cat(size,"")
      cat(mean(meanccorrRealZ),"")
      cat(mean(meanccorrPseudoZ),"")
      cat(sd(meanccorrRealZ),"")
      cat(sd(meanccorrPseudoZ),"")
      cat(length(meanccorrRealZ),"")
      cat(100*k1/(length(meanccorrRealZ)),"")
      cat(numberEpochen,"") #10
      # cat(numberEpochen*(numberEpochen-1),"")
      cat((mean(meanccorrRealZ)-mean(meanccorrPseudoZ))/sd(meanccorrPseudoZ),"")
      cat(mean(meanccorrRealZ[1:maxlagHz]),"")
      cat(mean(meanccorrRealZ[(maxlagHz+2):lagtimes2]),"") #14
      # cat(mean(meanccorrPseudoZ[1:maxlagHz]),"")
      # cat(mean(meanccorrPseudoZ[(maxlagHz+2):lagtimes2]),"") #16
      cat((mean(meanccorrRealZ[1:maxlagHz])-mean(meanccorrPseudoZ[1:maxlagHz]))/sd(meanccorrPseudoZ[1:maxlagHz]),"")
      cat((mean(meanccorrRealZ[(maxlagHz+2):lagtimes2])-mean(meanccorrPseudoZ[(maxlagHz+2):lagtimes2]))/sd(meanccorrPseudoZ[(maxlagHz+2):lagtimes2]),"") #18
      if (corr.no.abs) {
        cat(mean(meanccorrRealZNotAbs[which(meanccorrRealZNotAbs >= 0)]),"")
        cat(mean(meanccorrRealZNotAbs[which(meanccorrRealZNotAbs < 0)]),"")
        cat(length(meanccorrRealZNotAbs[which(meanccorrRealZNotAbs >= 0)]),"")
        cat(length(meanccorrRealZNotAbs[which(meanccorrRealZNotAbs < 0)]),"")
        cat(mean(meanccorrRealZNotAbs),"") #25
        cat(mean(meanccorrPseudoZNotAbs),"")
        cat(100*k1NotAbs/(length(meanccorrRealZNotAbs)),"")
        cat((mean(meanccorrRealZNotAbs)-mean(meanccorrPseudoZNotAbs))/sd(meanccorrPseudoZNotAbs),"")
      }
      cat("\n")
    }
    lapply(x, print1)
  }
  invisible(x)
}
