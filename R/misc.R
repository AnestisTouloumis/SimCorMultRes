check_cluster_size <- function(clsize){
  if (all.equal(clsize, as.integer(clsize)) != TRUE | clsize < 2)
    stop("'clsize' must be a positive integer greater than or equal to two")
}

check_ncategories <- function(ncategories){
  if (!is.numeric(ncategories) | ncategories < 3)
    stop("'ncategories' must be greater than or equal to three")
  ncategories <- as.integer(ncategories)
  if (all.equal(ncategories, as.integer(ncategories)) != TRUE | ncategories < 3)
    stop("'ncategories' must be a positive integer greater than or equal to three")
  ncategories
}

check_correlation_matrix <- function(cor.matrix, clsize, ncategories = 3,
                                     rfctn="rbin"){
  if (!is.numeric(cor.matrix))
    stop("'cor.matrix' must be numeric")
  if (!is.matrix(cor.matrix))
    stop("'cor.matrix' must be a matrix")
  if(rfctn=="rbin" | rfctn=="rmult.clm"){
    if (ncol(cor.matrix) != clsize | nrow(cor.matrix) != clsize)
      stop("'cor.matrix' must be a ", clsize, "x", clsize, " matrix")
  } else {
    ncateg2 <- ifelse(rfctn=="rmult.bcl", ncategories, ncategories-1)
    dimcor <- clsize * ncateg2
    if (ncol(cor.matrix) != dimcor | nrow(cor.matrix) != dimcor)
      stop("'cor.matrix' must be a ", dimcor, "x", dimcor, " matrix")
    for (i in 1:clsize) {
      diag.index <- 1:ncateg2 + (i - 1) * ncateg2
      cor.matrix[diag.index, diag.index] <- diag(1, ncateg2)
    }
  }
  if (!isSymmetric(cor.matrix))
    stop("'cor.matrix' must be a symmetric matrix")
  if (any(diag(cor.matrix) != 1))
    stop("the diagonal elements of 'cor.matrix' must be one")
  if (any(cor.matrix > 1) | any(cor.matrix < -1))
    stop("all the elements of 'cor.matrix' must be on [-1,1]")
  if (any(eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <= 0))
    stop("'cor.matrix' must be positive definite")
  if(rfctn=="rbin" | rfctn=="rmult.clm"){
    if(any(eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop("'cor.matrix' must be positive definite")
    } else {
      if(any(eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <= 0))
        stop("'cor.matrix' must respect the local independence of the alternatives and must be positive definite")
      }
    cor.matrix
  }

check_xformula <- function(xformula){
  lpformula <- as.formula(xformula)
  if (length(paste0(attr(terms(lpformula), "variables"))) == 1)
    stop("No covariates were found in 'formula' ")
  if (attr(terms(lpformula), "intercept") == 0) {
    lpformula <- update(lpformula, ~. + 1)
  }
  lpformula
}

check_intercepts <- function(intercepts, clsize, R=clsize, rfctn="rbin"){
  if (!is.numeric(intercepts))
    stop("'intercepts' must be numeric")
  if(rfctn=="rbin"){
  if (!(is.vector(intercepts) & !is.list(intercepts)))
    stop("'intercepts' must be a vector")
  if (length(intercepts) == 1)
    intercepts <- rep(intercepts, clsize)
  if (length(intercepts) != clsize)
    stop("'intercepts' must have either one or ", clsize, " elements")
  intercepts <- cbind(-Inf, intercepts, Inf)
  } else{
    if (
      !(is.vector(intercepts) & !is.list(intercepts)) & !is.matrix(intercepts))
      stop("'intercepts' must be a vector or a matrix")
    if (is.vector(intercepts) & !is.list(intercepts)) {
      if (length(intercepts) == 1)
        stop("'intercepts' must have at least 2 elements")
      if(rfctn=="rmult.clm" & any(diff(intercepts) <= 0))
        stop("'intercepts' must be increasing")
      ncategories <- length(intercepts) + 1
      if(rfctn=="rmult.clm"){
        intercepts <- matrix(intercepts, clsize, ncategories - 1, TRUE)
        intercepts <- cbind(-Inf, intercepts, Inf)
      } else {
        intercepts <- matrix(intercepts, R, clsize * (ncategories - 1),
                             byrow = TRUE)
        }
    } else {
      ncategories <- ncol(intercepts) + 1
      if(rfctn=="rmult.clm"){
      intercepts <- cbind(-Inf, intercepts, Inf)
      for (i in 1:clsize) {
        if (any(diff(intercepts[i, ]) <= 0))
          stop("'intercepts' must be increasing at each row")
      }
      } else {
        ncategories <- ncol(intercepts) + 1
        intercepts <- matrix(t(intercepts), R, clsize * (ncategories -
                                                           1), byrow = TRUE)
    }
    }
  }
  intercepts
}

check_betas <- function(betas, clsize)
{
  if (!(is.vector(betas) & !is.list(betas)) & !is.matrix(betas))
    stop("'betas' must be a vector or a matrix")
  if (!is.numeric(betas))
    stop("'betas' must be numeric")
  if (is.vector(betas) & !is.list(betas)) {
    betas <- rep(betas, clsize)
  } else {
  if (nrow(betas) != clsize)
        stop("The number of rows in 'betas' should be equal to 'clsize'")
      betas <- c(t(betas))
    }
  betas
}


create_linear_predictor <- function(betas, clsize, lpformula, xdata, ncategories=3,
                                    rfctn="rbin")
{
  Xmat <- model.matrix(lpformula, data = xdata)
  if(rfctn=="rmult.bcl"){
    Xmat <- apply(Xmat, 2, function(x) rep(x, each = ncategories))
    if (length(betas) != (clsize * ncategories * ncol(Xmat)))
      stop("The length of 'betas' does not match with the provided covariates")
  } else {
    Xmat <- matrix(Xmat[, -1], ncol = ncol(Xmat) - 1)
    if (length(betas) != (clsize) * ncol(Xmat))
      stop("The length of 'betas' does not match with the number of covariates")
  }
  lin.pred <- matrix(betas, nrow = nrow(Xmat), ncol = ncol(Xmat),
                     byrow = TRUE) * Xmat
  if(rfctn=="rmult.bcl"){
    lin.pred <- matrix(rowSums(lin.pred), ncol = ncategories * clsize,
                       byrow = TRUE)
  } else {
    lin.pred <- matrix(rowSums(lin.pred), ncol = clsize, byrow = TRUE)
  }
  as.matrix(lin.pred)
  }

 create_rlatent <- function(rlatent, R, link, clsize, ncategories = 3, rfctn="rbin"){
   if (is.null(rlatent)) {
     if (length(link) != 1)
       stop("The length of 'link' must be one")
       links <- c("probit", "logit", "cloglog", "cauchit")
       if (!is.element(link, links))
         stop("'link' must be 'probit','logit','cloglog' and/or 'cauchit'")
       distr <- switch(link, probit = "qnorm", logit = "qlogis",
                       cloglog = "qgumbel", cauchit = "qcauchy")
   cor.matrix <- check_correlation_matrix(cor.matrix, clsize, ncategories, rfctn)
   err <- rnorta(R, cor.matrix, rep(distr, nrow(cor.matrix)))
   if (distr == "qgumbel" & rfctn!="rmult.bcl") err <- -err
   } else {
     if (!is.matrix(rlatent))
       stop("'rlatent' must be a matrix")
     if (!is.numeric(rlatent))
       stop("'rlatent' must be numeric")
     if(rfctn=="rbin" | rfctn=="rmult.clm"){
       ncol_rlatent <- clsize
     } else if(rfctn=="rmult.bcl"){
       ncol_rlatent <- clsize * ncategories
     } else {
       ncol_rlatent <- clsize * (ncategories - 1)
     }
     ncol_rlatent
     if (nrow(rlatent) != R | ncol(rlatent) != ncol_rlatent)
       stop("'rlatent' must be a ", R, "x", ncol_rlatent, " matrix")
     cor.matrix <- NULL
     err <- rlatent
   }
   err
 }
