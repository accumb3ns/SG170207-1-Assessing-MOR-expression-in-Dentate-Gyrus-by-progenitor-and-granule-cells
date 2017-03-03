KLD2D_debug <- function (edge, bgIndex, fgIndex, expData, method = c("integers", 
                                                      "bandwidth")[2]) 
{
  require(KernSmooth)
  edge<<-edge
  xlevels <-expData[edge[1], c(bgIndex, fgIndex)]
  ylevels <-expData[edge[2], c(bgIndex, fgIndex)]

  x <- rank(expData[edge[1], c(bgIndex, fgIndex)], ties.method = "random")
  y <- rank(expData[edge[2], c(bgIndex, fgIndex)], ties.method = "random")
  N <- length(x)
  bgI <- 1:length(bgIndex)
  fgI <- length(bgIndex) + (1:length(fgIndex))
  
  ##if(!bw.nrd(x[fgI]) || !bw.nrd(y[fgI]) || !bw.nrd(x[bgI]) || !bw.nrd(x[bgI])) cat("skipped ", runs, " where\n", expData[edge[1], c(bgIndex, fgIndex)], "and \n", expData[edge[2], c(bgIndex, fgIndex)], "\n")
  ##else{
  fgWidth <<- c(bw.nrd(x[fgI]), bw.nrd(y[fgI]))
  bgWidth <<- c(bw.nrd(x[bgI]), bw.nrd(x[bgI]))
  gridSize <<- switch(method, integers = c(N, N), bandwidth = ceiling(N/c(min(fgWidth[1], 
                                                                              bgWidth[1]), min(fgWidth[2], bgWidth[2]))))
  ranges <- list(x = c(1, N), y = c(1, N))
  fgSmooth <- bkde2D(x = cbind(x[fgI], y[fgI]), bandwidth = fgWidth, 
                     range.x = ranges, gridsize = gridSize)
  fgP <- fgSmooth$fhat
  bgSmooth <- bkde2D(x = cbind(x[bgI], y[bgI]), bandwidth = bgWidth, 
                     range.x = ranges, gridsize = gridSize)
  bgP <- bgSmooth$fhat
  fgP <- pmax(fgP, 1e-20)
  bgP <- pmax(bgP, 1e-20)
  fgP <- fgP/sum(fgP)
  bgP <- bgP/sum(bgP)
  return((sum(fgP * log(fgP/bgP)) + sum(bgP * log(bgP/fgP)))/2)
}


runDemand_debug <- function (x, fgIndex = NULL, bgIndex = NULL, verbose = TRUE, 
          method = "bandwidth", keepLeaves = FALSE, alpha = 0.05) 
{
  if (is.null(fgIndex) | is.null(bgIndex)) 
    stop("Please provide sample (column of expression data) indices of case/control samples")
  if (any(is.na(c(fgIndex, bgIndex)))) {
    warning("Case indices contain NA values. These values are ignored")
    fgIndex <- as.vector(na.exclude(fgIndex))
    bgIndex <- as.vector(na.exclude(bgIndex))
  }
  if (length(fgIndex) < 3 | length(bgIndex) < 3) 
    stop("The number of samples in each class should be at least three")
  if (length(fgIndex) < 6 | length(bgIndex) < 6) 
    warning("DeMAND requires six samples (in each class) for optimal performance")
  expData <- x@exp
  if (any(is.na(expData))) 
    stop("Expression data contains NA values")
  if (any(is.infinite(x@exp))) 
    warning("Expression data contains infinite values")
  annot <- x@anno[, 2]
  if (any(is.na(annot))) 
    warning("Annotation data contains NA values")
  inputNetwork <- x@network
  if (any(is.na(inputNetwork[, 1:2]))) {
    warning("The network contains NA values, removing those lines")
    inputNetwork <- inputNetwork[apply(inputNetwork, 1, function(x) !any(is.na(x))), 
                                 ]
    x@network <- inputNetwork
  }
  platformGene <- unique(annot)
  vmsg <- function(x, verb = verbose) {
    if (verb) 
      message(x)
  }
  vmsg("Pruning the network")
  edgesToKeep <- apply(inputNetwork[, 1:2], 1, function(gg) all(gg %in% 
                                                                  platformGene))
  interactome <- inputNetwork[edgesToKeep, ]
  rm(edgesToKeep)
  analGene <- intersect(unique(as.vector(interactome[, 1:2])), 
                        platformGene)
  interactome[, 1:2] <- t(apply(interactome[, 1:2], 1, function(gg) c(min(gg), 
                                                                      max(gg))))
  dups <- duplicated(interactome)
  interactome <- interactome[!dups, 1:2]
  ppi <- if ("ppi" %in% colnames(inputNetwork)) {
    as.numeric(inputNetwork[!dups, "ppi"]) == 1
  }
  else {
    rep(F, sum(!dups))
  }
  vmsg("Keeping best probe per gene")
  CV <- sqrt(rowMeans(expData^2) - rowMeans(expData)^2)/rowMeans(expData)
  oGenes <- order(annot, -CV)
  dups <- duplicated(annot[oGenes])
  expData <- expData[oGenes[!dups], ]
  row.names(expData) <- annot[oGenes[!dups]]
  expData <- expData[analGene, ]
  getKLDpvalue <- function(kld, nullKLD) {
    rs <- sum(nullKLD >= kld)/length(nullKLD)
    return(min(1, rs))
  }
  vmsg("Make a null distribution for KL divergence.....")
  p1 <- sample(analGene, min(max(length(analGene), 1000), 10000), 
               replace = T)
  p2 <- sample(analGene, min(max(length(analGene), 1000), 10000), 
               replace = T)
  pKeep <- !(p1 == p2)
  permuteInteractome <- cbind(p1[pKeep], p2[pKeep])
  permuteInteractome_initialized <<- permuteInteractome
  permuteInteractome <- t(apply(permuteInteractome, 1, function(x) c(min(x), 
                                                                     max(x))))
  permuteInteractome_postMinMax <<- permuteInteractome
  dups <- duplicated(permuteInteractome)
  nullBgIndex <- bgIndex
  nullFgIndex <- fgIndex
  nullKLD <- apply(permuteInteractome, 1, KLD2D_debug, nullBgIndex, 
                   nullFgIndex, expData, method)
  vmsg("Measure dysregulation of the interactions.....")
  KLDmat <- apply(interactome, 1, KLD2D_debug, bgIndex, fgIndex, 
                  expData, method)
  pfit <- pareto.fit_debug(data = nullKLD, threshold = quantile(nullKLD, 
                                                          probs = 1 - alpha))
  KLDpvec <- ppareto(x = KLDmat, threshold = pfit$xmin, exponent = pfit$exponent, 
                     lower.tail = F) * alpha
  pToReplace <- KLDmat < quantile(nullKLD, probs = 1 - alpha)
  KLDpvec[pToReplace] <- sapply(X = KLDmat[pToReplace], FUN = function(K) getKLDpvalue(K, 
                                                                                       nullKLD))
  edgeKLD <- cbind(interactome, KLDmat, KLDpvec)
  colnames(edgeKLD) <- c("gene1", "gene2", "KLD", "KLD.p")
  vmsg("Estimate dysregulation of the genes.....")
  intPval <- apply(as.matrix(analGene), 1, integratePvalues, 
                   edgeKLD, expData, ppi, keepLeaves)
  intPvalAdjustedp <- p.adjust(intPval, "fdr")
  finalPrediction <- data.frame(moaGene = analGene, Pvalue = intPval, 
                                FDR = intPvalAdjustedp)
  finalPrediction <- finalPrediction[order(intPval, decreasing = F), 
                                     ]
  rownames(finalPrediction) <- 1:dim(finalPrediction)[1]
  x@moa <- finalPrediction
  x@KLD <- as.data.frame(edgeKLD[order(KLDpvec), ])
  return(x)
}

pareto.fit_debug <- function (data, threshold, method = "ml") 
{
  if (threshold == "find") {
    return(.pareto.fit.threshold(data, method = method))
  }
  switch(method, ml = {
    return(.pareto.fit.ml_debug(data, threshold))
  }, regression.cdf = {
    return(.pareto.fit.regression.cdf(data, threshold))
  }, {
    cat("Unknown method\n")
    return(NA)
  })
}

.pareto.fit.ml_debug <- function (data, threshold) 
{
  data <- data[data >= threshold]
  n <- length(data)
  x <- data/threshold
  alpha <- 1 + n/sum(log(x))
  loglike = pareto.loglike_debug(data, threshold, alpha)
  ks.dist <- .ks.dist.fixed.pareto_debug(data, threshold = threshold, 
                                   exponent = alpha)
  fit <- list(type = "pareto", exponent = alpha, xmin = threshold, 
              loglike = loglike, ks.dist = ks.dist, samples.over.threshold = n)
  return(fit)
}

pareto.loglike_debug <- function (x, threshold, exponent) 
{
  L <- sum(dpareto_debug(x, threshold = threshold, exponent = exponent, 
                   log = TRUE))
  return(L)
}

.ks.dist.fixed.pareto_debug <- function (data, threshold, exponent) 
{
  data <- data[data >= threshold]
  d <- suppressWarnings(ks.test(data, ppareto, threshold = threshold, 
                                exponent = exponent))
  return(as.vector(d$statistic))
}

dpareto_debug <- function (x, threshold = 1, exponent, log = FALSE) 
{
  if (!log) {
    prefactor <- (exponent - 1)/threshold
    f <- function(x) {
      prefactor * (x/threshold)^(-exponent)
    }
  }
  else {
    prefactor.log <- log(exponent - 1) - log(threshold)
    f <- function(x) {
      prefactor.log - exponent * (log(x) - log(threshold))
    }
  }
  d <- ifelse(x < threshold, NA, f(x))
  return(d)
}

ppareto <- function (x, threshold = 1, exponent, lower.tail = TRUE, log.p = FALSE) 
{
  if ((!lower.tail) && (!log.p)) {
    f <- function(x) {
      (x/threshold)^(1 - exponent)
    }
  }
  if ((lower.tail) && (!log.p)) {
    f <- function(x) {
      1 - (x/threshold)^(1 - exponent)
    }
  }
  if ((!lower.tail) && (log.p)) {
    f <- function(x) {
      (1 - exponent) * (log(x) - log(threshold))
    }
  }
  if ((lower.tail) && (log.p)) {
    f <- function(x) {
      log(1 - (x/threshold)^(1 - exponent))
    }
  }
  p <- ifelse(x < threshold, NA, f(x))
  return(p)
}

integratePvalues <- function (g, network, expData, ppi, keepLeaves) 
{
  neighborEdges <- which(network[, 1] %in% g | network[, 2] %in% 
                           g)
  N <- length(neighborEdges)
  if (N < 2) 
    if (keepLeaves) {
      return(as.numeric(network[neighborEdges, 4][1]))
    }
  else {
    return(1)
  }
  pvals <- as.numeric(as.vector(network[neighborEdges, 4]))
  pvals <- pmax(pvals, 1e-20)
  fisherChisq <- -2 * sum(log(pvals))
  neighborGenes <- apply(network[neighborEdges, 1:2], 1, function(x) setdiff(x, 
                                                                             g))
  if (length(unique(neighborGenes)) < 2) 
    if (keepLeaves) {
      return(as.numeric(network[neighborEdges, 4][1]))
    }
  else {
    return(1)
  }
  neighborExp <- expData[neighborGenes, ]
  gExp <- expData[g, ]
  resids <- lm(t(neighborExp) ~ gExp)$residuals
  resids[, ppi[neighborEdges]] <- t(expData[neighborGenes, 
                                            ][ppi[neighborEdges], ])
  covMat <- cov(resids)
  covMat <- covMat[lower.tri(covMat)]
  I <- covMat < 0
  covMat[I] <- covMat[I] * (3.27 + 0.71 * covMat[I])
  covMat[!I] <- covMat[!I] * (3.25 + 0.75 * covMat[!I])
  c <- 4 * N/(4 * N + 2 * sum(covMat))
  dFreedom <- 2 * N * c
  correctedChisq <- fisherChisq * c
  return(pchisq(q = correctedChisq, df = dFreedom, lower.tail = F))
}
