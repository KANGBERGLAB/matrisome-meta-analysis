function(Exp1, Pseudotime1, PseudotimeRange1 = NULL,
                        SmoothLength1 = 100, ByBin1 = c("Equal.Pseudotime", "Equal.Cells")) {
  if (ncol(Exp1) != nrow(Pseudotime1)) {
    stop("The length of pseudotime is not equal to the number of cells")
  } else if (!is.element("Pseudotime", colnames(Pseudotime1))) {
    stop("No pseudotime inforamtion in variable Pseudotime1")
  }

  if (ByBin1[1] == "Equal.Pseudotime") {
    if (is.null(PseudotimeRange1)) {
      PseudotimeBin1 <- seq(min(Pseudotime1$Pseudotime), max(Pseudotime1$Pseudotime),
                            length.out = SmoothLength1 + 1)
    } else {
      PseudotimeBin1 <- seq(PseudotimeRange1[1], PseudotimeRange1[2],
                            length.out = SmoothLength1 + 1)
    }
  } else {
    Pseudotime1 <- Pseudotime1[order(Pseudotime1$Pseudotime), ]
    Bin1 <- ceiling(nrow(Pseudotime1) / SmoothLength1)
    PseudotimeBin1 <- seq(1, nrow(Pseudotime1), by = Bin1)
    PseudotimeBin1[length(PseudotimeBin1) + 1] <- nrow(Pseudotime1)
  }

  Exp2 <- array(0, dim = c(nrow(Exp1), SmoothLength1))
  for (i in 1:(length(PseudotimeBin1) - 1)) {
    if (ByBin1[1] == "Equal.Pseudotime") {
      if (i == length(PseudotimeBin1) - 1) {
        Cells1 <- Pseudotime1[Pseudotime1$Pseudotime >= PseudotimeBin1[i] &
                                Pseudotime1$Pseudotime <= PseudotimeBin1[i + 1], ]
      } else {
        Cells1 <- Pseudotime1[Pseudotime1$Pseudotime >= PseudotimeBin1[i] &
                                Pseudotime1$Pseudotime < PseudotimeBin1[i + 1], ]
      }
    } else {
      if (i == length(PseudotimeBin1) - 1) {
        Cells1 <- Pseudotime1[PseudotimeBin1[i]:PseudotimeBin1[i + 1], ]
      } else {
        Cells1 <- Pseudotime1[PseudotimeBin1[i]:(PseudotimeBin1[i + 1] - 1), ]
      }
    }

    if (nrow(Cells1) > 1) {
      Exp2[, i] <- rowMeans(Exp1[, match(rownames(Cells1), colnames(Exp1))])
    } else if (nrow(Cells1) == 1) {
      Exp2[, i] <- Exp1[, match(rownames(Cells1), colnames(Exp1))]
    }
  }
  rownames(Exp2) <- rownames(Exp1)

  return(Exp2)
}
