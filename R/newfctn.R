rbin2 <- function(clsize = clsize, intercepts = intercepts, betas = betas,
                 xformula = formula(xdata), xdata = parent.frame(), link = "logit",
                 cor.matrix = cor.matrix, rlatent = NULL) {
  check_cluster_size(clsize)
  intercepts <- check_intercepts(intercepts, clsize, clsize, "rbin")
  betas <- check_betas(betas, clsize)
  lpformula <- check_xformula(xformula)
  lin.pred <- create_linear_predictor(betas, clsize, lpformula, xdata)
  R <- nrow(lin.pred)
  rlatent <- create_rlatent(rlatent, R, link, clsize)
  U <- lin.pred + rlatent
  Ysim <- matrix(0, R, clsize)
  for (i in 1:clsize) Ysim[, i] <-
    cut(U[, i] - 2 * lin.pred[, i], intercepts[i, ], labels = FALSE)
  Ysim <- 2 - Ysim
  y <- c(t(Ysim))
  id <- rep(1:R, each = clsize)
  time <- rep(1:clsize, R)
  rownames(Ysim) <- rownames(rlatent) <- paste("i", 1:R, sep = "=")
  colnames(Ysim) <- colnames(rlatent) <- paste("t", 1:clsize, sep = "=")
  simdata <- data.frame(y, model.frame(formula = lpformula, data = xdata), id, time)
  list(Ysim = Ysim, simdata = simdata, rlatent = rlatent)
}


  rmult.bcl2 <- function(clsize = clsize, ncategories = ncategories,
                      betas = betas, xformula = formula(xdata),
                      xdata = parent.frame(), cor.matrix = cor.matrix,
                      rlatent = NULL) {
  check_cluster_size(clsize)
  ncategories <- check_ncategories(ncategories)
  betas <- check_betas(betas, clsize)
  lpformula <- check_xformula(xformula)
  lin.pred <- create_linear_predictor(betas,clsize,lpformula,xdata,ncategories,"rmult.bcl")
  R <- nrow(lin.pred)
  rlatent <- create_rlatent(rlatent,R,"cloglog", clsize, ncategories,"rmult.bcl")
  U <- lin.pred + rlatent
  U <- matrix(as.vector(t(U)), nrow = clsize * R, ncol = ncategories,
              TRUE)
  Ysim <- apply(U, 1, which.max)
  Ysim <- matrix(Ysim, ncol = clsize, byrow = TRUE)
  id <- rep(1:R, each = clsize)
  time <- rep(1:clsize, R)
  y <- c(t(Ysim))
  lpformula <- update(lpformula, ~. - 1)
  rownames(Ysim) <- rownames(rlatent) <- paste("i", 1:R, sep = "=")
  colnames(Ysim) <- paste("t", 1:clsize, sep = "=")
  colnames(rlatent) <- paste("t=", rep(1:clsize, each = ncategories), " & j=",
                         rep(1:ncategories, clsize), sep = "")
  simdata <- data.frame(y, model.frame(formula = lpformula, data = xdata), id, time)
  list(Ysim = Ysim, simdata = simdata, rlatent = rlatent)
}


  rmult.clm2 <- function(clsize = clsize, intercepts = intercepts, betas = betas,
                        xformula = formula(xdata), xdata = parent.frame(), link = "logit",
                        cor.matrix = cor.matrix, rlatent = NULL) {
    check_cluster_size(clsize)
    intercepts <- check_intercepts(intercepts,clsize,clsize,"rmult.clm")
    betas <- check_betas(betas,clsize)
    lpformula <- check_xformula(xformula)
    lin.pred <- create_linear_predictor(betas, clsize, lpformula, xdata, 3, "rmult.clm")
    R <- nrow(lin.pred)
    rlatent <- create_rlatent(rlatent,R, link, clsize, 3, "rmult.clm")
    U <- -lin.pred + rlatent
    Ysim <- matrix(0, R, clsize)
    for (i in 1:clsize) Ysim[, i] <-
      cut(U[, i], intercepts[i, ], labels = FALSE)
    id <- rep(1:R, each = clsize)
    time <- rep(1:clsize, R)
    y <- c(t(Ysim))
    rownames(Ysim) <- rownames(rlatent) <- paste("i", 1:R, sep = "=")
    colnames(Ysim) <- colnames(rlatent) <- paste("t", 1:clsize, sep = "=")
    simdata <- data.frame(y, model.frame(formula = lpformula, data = xdata),
                          id, time)
    list(Ysim = Ysim, simdata = simdata, rlatent = rlatent)
  }




  rmult.crm2 <- function(clsize = clsize, intercepts = intercepts, betas = betas,
                        xformula = formula(xdata), xdata = parent.frame(), link = "logit",
                        cor.matrix = cor.matrix, rlatent = NULL) {
    check_cluster_size(clsize)
    betas <- check_betas(betas,clsize)
    lpformula <- check_xformula(xformula)
    lin.pred <- create_linear_predictor(betas, clsize, lpformula, xdata, 3, "rmult.clm")
    R <- nrow(lin.pred)
    intercepts <- check_intercepts(intercepts,clsize,R,"rmult.crm")
    ncategories <- ncol(intercepts)/clsize + 1
    lin.pred.extended <- t(apply(lin.pred, 1,
                                 function(x) rep(x, each = ncategories - 1)))
    rlatent <- create_rlatent(rlatent,R,link,clsize,ncategories,"rmult.crm")
    U <- rlatent - lin.pred.extended
    Ysim <- matrix(as.numeric(t(U <= intercepts)), R * clsize, ncategories -
                     1, TRUE)
    for (i in 1:(ncategories - 1)) Ysim[, i] <- ifelse(Ysim[, i] == 1,
                                                       i, ncategories)
    Ysim <- apply(Ysim, 1, min)
    Ysim <- matrix(Ysim, R, clsize, byrow = TRUE)
    id <- rep(1:R, each = clsize)
    time <- rep(1:clsize, R)
    y <- c(t(Ysim))
    rownames(Ysim) <- rownames(rlatent) <- paste("i", 1:R, sep = "=")
    colnames(Ysim) <- paste("t", 1:clsize, sep = "=")
    colnames(rlatent) <- paste("t=", rep(1:clsize, each = ncategories - 1),
                           " & j=", rep(1:(ncategories - 1), clsize), sep = "")
    simdata <- data.frame(y, model.frame(formula = lpformula, data = xdata),
                          id, time)
    list(Ysim = Ysim, simdata = simdata, rlatent = rlatent)
  }
