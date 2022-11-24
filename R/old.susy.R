old.susy = function(
  filename="MEA-topstage-C3-Beethoven.txt", separator=" ", missings=".", datahead=TRUE,
  spalte1=2, spalte2=5, epoche=30, fps=15, maxlag=3, automatic=FALSE, plotid=3, plotid2=5,
  pseudosVereinfachen=FALSE, anzahlPseudosGesamt=500, zufallsdaten=FALSE, anzahlzufallsdaten=10000,
  korrOhneBetragAnzeigen=TRUE, permutation=FALSE
  ) {
  # SUSY surrogate analysis of dyadic synchrony
  # David Leander Tschacher, 2016
  # 2018: wenn plotid = 1 gewahlt wird, wird die CCF als file exportiert!
  # 2019: laengste ZR muss nicht in erste Spalte, missings stehen als Punkt "." im Text-file
  # 2020: Permutation der ZR

  if (!zufallsdaten) {
    data<-read.csv(file=filename,header=datahead,sep=separator,na.strings=missings)

    if (permutation) {
      if (!requireNamespace("gtools", quietly=TRUE))
        stop("old.susy, when provided permutation=TRUE, requires gtools, install gtools package and retry")
      list_combination<-gtools::combinations(n=ncol(data),r=2,v=colnames(data),set=TRUE,repeats.allowed=FALSE)
      comb<-rbind(list_combination[,1],list_combination[,2])
      data<-data[match(comb,colnames(data))]
    }
  }

  if (zufallsdaten) {
    automatic = FALSE
    a = runif(anzahlzufallsdaten,300,330)
    b = runif(anzahlzufallsdaten,300,330)
    variablenname1 = "randomA"
    variablenname2 = "randomB"
  }

  startColumn = 1
  numberColumns = 2

  if (automatic) {
    numberColumns = length(data[1,]) + 1 - startColumn
  }
  firstloop = TRUE

  while (numberColumns >= 2) {
    if (automatic) {
      a = data[, startColumn]
      b = data[, startColumn + 1]
      variablenname1 = colnames(data)[startColumn]
      variablenname2 = colnames(data)[startColumn + 1]
    } else {
      if (!zufallsdaten) {
        a = data[, spalte1]
        b = data[, spalte2]
        variablenname1 = colnames(data)[spalte1]
        variablenname2 = colnames(data)[spalte2]
      }
      numberColumns = 0
    }
    a = a[!is.na(a)]
    b = b[!is.na(b)]
    maxlag = maxlag*fps
    lagtimes2 = maxlag*2 + 1

    size = length(a)
    size2 = length(b)
    if (size2 < size) {
      size = size2
    }

    range = epoche * fps

    numberEpochen = round(size/range-0.499999)

    if (pseudosVereinfachen) {
      anzahlPseudosProEpoche = floor(anzahlPseudosGesamt/numberEpochen)
      if (anzahlPseudosProEpoche < 1){
        cat("error: anzahlPseudosProGesamt < numberEpochen, automatisch auf min vergrossert\n")
        anzahlPseudosProEpoche = 1
      }
      if (anzahlPseudosProEpoche >= numberEpochen){
        cat("error: anzahlPseudosGesamt zu gross, automatisch auf max verkleinert\n")
        anzahlPseudosProEpoche = numberEpochen-1
      }
      if (anzahlPseudosProEpoche == numberEpochen-1) {
        pseudosVereinfachen = FALSE
      }
    } else {
      anzahlPseudosProEpoche = numberEpochen-1
    }

    #epochen generieren
    eps = array(0:0, dim=c(2, numberEpochen, range))

    i0 = 1
    i2 = 1
    while (i0 <= numberEpochen) {
      i1 = 1
      while (i1 <= range) {
        eps[1,i0,i1] = a[i2]
        eps[2,i0,i1] = b[i2]
        i1 = i1 + 1
        i2 = i2 + 1
      }
      i0 = i0 + 1
    }

    #durchschnitte
    meansReal = seq(from = 0, to = 0, length.out = round(size/range-0.5))
    meansPseudo = seq(from = 0, to = 0, length.out = round(size/range-0.5))

    #meanccorr mit nullen fullen, lange = lagtimes2
    meanccorrReal = seq(from = 0, to = 0, length.out = lagtimes2)
    meanccorrPseudo = seq(from = 0, to = 0, length.out = lagtimes2)
    meanccorrRealZ = seq(from = 0, to = 0, length.out = lagtimes2)
    meanccorrRealZNotAbs = seq(from = 0, to = 0, length.out = lagtimes2)
    meanccorrPseudoZ = seq(from = 0, to = 0, length.out = lagtimes2)
    meanccorrPseudoZNotAbs = seq(from = 0, to = 0, length.out = lagtimes2)

    #meanccorr mit nullen fullen, lange = lagtimes2
    nReal = seq(from = 0, to = 0, length.out = numberEpochen)
    nPseudo = seq(from = 0, to = 0, length.out = numberEpochen)

    i = 1

    if (pseudosVereinfachen) {
      while (i <= numberEpochen) {
        result = ccf(eps[1,i,], eps[2,i,], lag.max = maxlag, plot = FALSE)
        if (!is.nan(result$acf[1])) {
          j = 1
          while (j <= lagtimes2) {
            meanccorrRealZ[j] = meanccorrRealZ[j] +
              abs(0.5*log((1+result$acf[j])/(1-result$acf[j])))
            meanccorrRealZNotAbs[j] = meanccorrRealZNotAbs[j] +
              0.5*log((1+result$acf[j])/(1-result$acf[j]))
            meanccorrReal[j] = meanccorrReal[j] + abs(result$acf[j])
            nReal[i] = nReal[i] + abs(result$acf[j])
            j = j + 1
          }
          nReal[i] = nReal[i]/lagtimes2
        }
        i = i + 1
      }

      isOcc = seq(from = 0, to = 0, length.out = numberEpochen)
      i = 1
      n = 0
      v = 1
      while (i <= numberEpochen) {
        h = 1
        while (h <= anzahlPseudosProEpoche) {
          repeat {
            v = floor(runif(1,1,numberEpochen+1))
            if (isOcc[v] != i && v != i){
              isOcc[v] = i
              break
            }
          }
          n = n + 1
          result = ccf(eps[1,i,], eps[2,v,], lag.max = maxlag, plot = FALSE)
          if (!is.nan(result$acf[1])) {
            j = 1
            while (j <= lagtimes2) {
              meanccorrPseudoZ[j] = meanccorrPseudoZ[j] +
                abs(0.5*log((1+result$acf[j])/(1-result$acf[j])))
              meanccorrPseudoZNotAbs[j] = meanccorrPseudoZNotAbs[j] +
                0.5*log((1+result$acf[j])/(1-result$acf[j]))
              meanccorrPseudo[j] = meanccorrPseudo[j] + abs(result$acf[j])
              nPseudo[i] = nPseudo[i] + abs(result$acf[j])
              j = j + 1
            }
          }
          h = h + 1
        }
        nPseudo[i] = nPseudo[i]/(anzahlPseudosProEpoche*lagtimes2)
        i = i + 1
      }
      meanccorrPseudo = meanccorrPseudo/n
      meanccorrPseudoZ = meanccorrPseudoZ/n
      meanccorrPseudoZNotAbs = meanccorrPseudoZNotAbs/n
    } else {
      while (i <= numberEpochen) {
        h = 1
        while (h <= numberEpochen) {
          result = ccf(eps[1,i,], eps[2,h,], lag.max = maxlag, plot = FALSE)
          if (!is.nan(result$acf[1])) {
            if (i == h) {
              j = 1
              while (j <= lagtimes2) {
                meanccorrRealZ[j] = meanccorrRealZ[j] +
                  abs(0.5*log((1+result$acf[j])/(1-result$acf[j])))
                meanccorrRealZNotAbs[j] = meanccorrRealZNotAbs[j] +
                  0.5*log((1+result$acf[j])/(1-result$acf[j]))
                meanccorrReal[j] = meanccorrReal[j] + abs(result$acf[j])
                nReal[i] = nReal[i] + abs(result$acf[j])
                j = j + 1
              }
            } else {
              j = 1
              while (j <= lagtimes2) {
                meanccorrPseudoZ[j] = meanccorrPseudoZ[j] +
                  abs(0.5*log((1+result$acf[j])/(1-result$acf[j])))
                meanccorrPseudoZNotAbs[j] = meanccorrPseudoZNotAbs[j] +
                  0.5*log((1+result$acf[j])/(1-result$acf[j]))
                meanccorrPseudo[j] = meanccorrPseudo[j] + abs(result$acf[j])
                nPseudo[i] = nPseudo[i] + abs(result$acf[j])
                j = j + 1
              }
            }
          }
          h = h + 1
        }
        nReal[i] = nReal[i]/lagtimes2
        nPseudo[i] = nPseudo[i]/((numberEpochen-1)*lagtimes2)
        i = i + 1
      }
      meanccorrPseudo = meanccorrPseudo/((numberEpochen-1)*numberEpochen)
      meanccorrPseudoZ = meanccorrPseudoZ/((numberEpochen-1)*numberEpochen)
      meanccorrPseudoZNotAbs = meanccorrPseudoZNotAbs/((numberEpochen-1)*numberEpochen)
    }

    meanccorrReal = meanccorrReal/numberEpochen
    meanccorrRealZ = meanccorrRealZ/numberEpochen
    meanccorrRealZNotAbs = meanccorrRealZNotAbs/numberEpochen

    t = 1
    k1 = 0
    k1NotAbs = 0
    while (t <= lagtimes2) {
      if (meanccorrRealZ[t] > meanccorrPseudoZ[t]) {
        k1 = k1 + 1
      }
      if (meanccorrRealZNotAbs[t] > meanccorrPseudoZNotAbs[t]) {
        k1NotAbs = k1NotAbs + 1
      }
      t = t + 1
    }

    numberColumns = numberColumns - 2
    startColumn = startColumn + 2

    if (!automatic) {
      #zeichne plot
      zaehler3 = 0
      while (zaehler3 < 2) {
        if (plotid == 1) {

          #synchronie meanccorrReal meanccorrPseudo
          min0 = min(meanccorrReal, meanccorrPseudo)
          max0 = max(meanccorrReal, meanccorrPseudo)
          title0 =  paste("Synchrony",variablenname1,variablenname2," segment:",epoche,"s;",anzahlPseudosProEpoche*numberEpochen,"pseudos")
          plot(seq(from = -maxlag , to = maxlag), meanccorrReal, ylim=c(min0, max0+ 0.19*(max0-min0)), main=title0,
               xlab="lag", ylab="correlation", type="l", col="green", lwd = 4)
          lines(seq(from = -maxlag , to = maxlag), meanccorrPseudo, col="red", lwd = 4)
          legend("topright", pch = c(3), col = c("green", "red"),legend = c("meanccorr", "meanccorr pseudo"))

          dfrm <- data.frame(xdata = seq(from = -maxlag , to = maxlag), meanccorrPseudo = meanccorrPseudo, meanccorrReal = meanccorrReal)
          write.table(dfrm, file="CrossCorrelations.txt", sep=",", row.names=FALSE, col.names=TRUE)

        } else if (plotid == 2) {

          #Epochensynchronie nReal nPseudo
          title0 =  paste("Segmentsynchronie",variablenname1,variablenname2," segment:",epoche,"s;",anzahlPseudosProEpoche*numberEpochen,"pseudos")
          plot(seq(from = 1 , to = numberEpochen ), nReal, ylim=c(min(nReal, nPseudo), max(nReal, nPseudo)), main=title0,
               xlab="Segment-Nr.", ylab="correlation", type="l", col="green")
          lines(seq(from = 1, to = numberEpochen ), nPseudo, col="red", lwd = 2)
          legend("topright", pch = c(3), col = c("green", "red"),legend = c("cor real", "cor pseudo"))

        } else if (plotid == 3) {

          #Synchrony meanccorrRealZ meanccorrPseudoZ
          min0 = min(meanccorrRealZ, meanccorrPseudoZ)
          max0 = max(meanccorrRealZ, meanccorrPseudoZ)
          title0 =  paste("Z-Synchrony",variablenname1,variablenname2," segment:",epoche,"s;",anzahlPseudosProEpoche*numberEpochen,"pseudos")
          plot(seq(from = -maxlag , to = maxlag), meanccorrRealZ, ylim=c(min0, max0+ 0.19*(max0-min0)),
               main=title0, xlab="lag", ylab="correlation", type="l", col="green", lwd = 4)
          lines(seq(from = -maxlag , to = maxlag), meanccorrPseudoZ, col="red", lwd = 4)
          legend("topright", pch = c(3), col = c("green", "red"),legend = c("Z(meanccorr)", "Z(meanccorr pseudo)"))

        } else if (plotid == 4) {

          #Zeige die Verteilung der data an
          title0 = paste("eingelesene Zeitreihen")
          plot(seq(from = 1 , to = size), a,ylim=c(min(a, b), max(a, b)), main=title0,
               xlab="Zeit", ylab="Wert", type="l", pch=20, col="green")
          points(seq(from = 1 , to = size), b, type="l", pch=20, col="red", lwd = 2)
          legend("topright", pch = c(20,20), col = c("green", "red"),legend = c(colnames(data)[spalte1], 		colnames(data)[spalte2]))

        } else if (plotid == 5) {

          #Synchrony meanccorrRealZNotAbs meanccorrPseudoZNotAbs
          min0 = min(meanccorrRealZNotAbs, meanccorrPseudoZNotAbs)
          max0 = max(meanccorrRealZNotAbs, meanccorrPseudoZNotAbs)
          title0 =  paste("Z-Synchrony NotAbs",variablenname1,variablenname2," segment:",epoche,"s;",anzahlPseudosProEpoche*numberEpochen,"pseudos")
          plot(seq(from = -maxlag , to = maxlag), meanccorrRealZNotAbs, ylim=c(min0, max0+ 0.19*(max0-min0)),
               main=title0, xlab="lag", ylab="correlation", type="l", col="green", lwd = 4)
          lines(seq(from = -maxlag , to = maxlag), meanccorrPseudoZNotAbs, col="red", lwd = 4)
          legend("topright", pch = c(3), col = c("green", "red"),legend = c("Z(meanccorr NotAbs)", "Z(meanccorr pseudo NotAbs)"))

        }
        zaehler3 = zaehler3 + 1
        if (plotid2 == 0){
          zaehler3 = zaehler3 + 1
        } else {
          plotid = plotid2
          plotid2 = 0
          dev.new()
        }
      }
    }

    if (firstloop) {
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
      if (korrOhneBetragAnzeigen) {
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
      firstloop = FALSE
    }

    {
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
      cat(mean(meanccorrRealZ[1:maxlag]),"")
      cat(mean(meanccorrRealZ[(maxlag+2):lagtimes2]),"") #14
      # cat(mean(meanccorrPseudoZ[1:maxlag]),"")
      # cat(mean(meanccorrPseudoZ[(maxlag+2):lagtimes2]),"") #16
      cat((mean(meanccorrRealZ[1:maxlag])-mean(meanccorrPseudoZ[1:maxlag]))/sd(meanccorrPseudoZ[1:maxlag]),"")
      cat((mean(meanccorrRealZ[(maxlag+2):lagtimes2])-mean(meanccorrPseudoZ[(maxlag+2):lagtimes2]))/sd(meanccorrPseudoZ[(maxlag+2):lagtimes2]),"") #18
      if (korrOhneBetragAnzeigen) {
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
  }
}
