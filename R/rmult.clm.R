rmult.clm <- function(clsize = clsize, intercepts = intercepts, betas = betas, xformula = formula(xdata), 
    xdata = parent.frame(), link = "logit", cor.matrix = cor.matrix, rlatent = NULL) {
    if (all.equal(clsize, as.integer(clsize)) != TRUE | clsize < 2) 
        stop("'clsize' must be a positive integer greater than or equal to two")
    if (!(is.vector(intercepts) & !is.list(intercepts)) & !is.matrix(intercepts)) 
        stop("'intercepts' must be a vector or a matrix")
    if (!is.numeric(intercepts)) 
        stop("'intercepts' must be numeric")
    if (is.vector(intercepts) & !is.list(intercepts)) {
        if (length(intercepts) == 1) 
            stop("'intercepts' must have at least 2 elements")
        if (any(diff(intercepts) <= 0)) 
            stop("'intercepts' must be increasing")
        ncategories <- length(intercepts) + 1
        intercepts <- matrix(intercepts, clsize, ncategories - 1, TRUE)
        intercepts <- cbind(-Inf, intercepts, Inf)
    } else {
        ncategories <- ncol(intercepts) + 1
        intercepts <- cbind(-Inf, intercepts, Inf)
        for (i in 1:clsize) {
            if (any(diff(intercepts[i, ]) <= 0)) 
                stop("'intercepts' must be increasing at each row")
        }
    }
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
    lpformula <- as.formula(xformula)
    if (length(paste0(attr(terms(lpformula), "variables"))) == 1) 
        stop("No covariates were found in 'formula' ")
    Xmat <- model.matrix(lpformula, data = xdata)
    if (attr(terms(lpformula), "intercept") == 0) {
        lpformula <- update(lpformula, ~. + 1)
    }
    Xmat <- matrix(Xmat[, -1], ncol = ncol(Xmat) - 1)
    if (length(betas) != (clsize) * ncol(Xmat)) 
        stop("The length of 'betas' does not match with the number of covariates")
    lin.pred <- matrix(betas, nrow = nrow(Xmat), ncol = ncol(Xmat), byrow = TRUE) * 
        Xmat
    lin.pred <- matrix(rowSums(lin.pred), ncol = clsize, byrow = TRUE)
    lin.pred <- as.matrix(lin.pred)
    R <- nrow(lin.pred)
    if (is.null(rlatent)) {
        if (length(link) != 1) 
            stop("The length of 'link' must be one")
        links <- c("probit", "logit", "cloglog", "cauchit")
        if (!is.element(link, links)) 
            stop("'link' must be 'probit','logit','cloglog' and/or 'cauchit'")
        distr <- switch(link, probit = "qnorm", logit = "qlogis", cloglog = "qgumbel", 
            cauchit = "qcauchy")
        if (!is.numeric(cor.matrix)) 
            stop("'cor.matrix' must be numeric")
        if (!is.matrix(cor.matrix)) 
            stop("'cor.matrix' must be matrix")
        if (ncol(cor.matrix) != clsize | nrow(cor.matrix) != clsize) 
            stop("'cor.matrix' must be a ", clsize, "x", clsize, " matrix")
        if (!isSymmetric(cor.matrix)) 
            stop("'cor.matrix' must be a symmetric matrix")
        if (any(diag(cor.matrix) != 1)) 
            stop("the diagonal elements of 'cor.matrix' must be one")
        if (any(cor.matrix > 1) | any(cor.matrix < -1)) 
            stop("all the elements of 'cor.matrix' must be on [-1,1]")
        if (any(eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <= 
            0)) 
            stop("'cor.matrix' must be positive definite")
        if (length(distr) == 1) 
            err <- rnorta(R, cor.matrix, rep(distr, clsize))
        if (distr == "qgumbel") 
            err <- -err
    } else {
        if (!is.matrix(rlatent)) 
            stop("'rlatent' must be a matrix")
        if (!is.numeric(rlatent)) 
            stop("'rlatent' must be numeric")
        if (ncol(rlatent) != clsize | nrow(rlatent) != R) 
            stop("'rlatent' must be a ", R, "x", clsize, " matrix")
        cor.matrix <- NULL
        err <- rlatent
    }
    U <- -lin.pred + err
    Ysim <- matrix(0, R, clsize)
    for (i in 1:clsize) Ysim[, i] <- cut(U[, i], intercepts[i, ], labels = FALSE)
    id <- rep(1:R, each = clsize)
    time <- rep(1:clsize, R)
    y <- c(t(Ysim))
    rownames(Ysim) <- rownames(err) <- paste("i", 1:R, sep = "=")
    colnames(Ysim) <- colnames(err) <- paste("t", 1:clsize, sep = "=")
    simdata <- data.frame(y, model.frame(formula = lpformula, data = xdata), id, 
        time)
    list(Ysim = Ysim, simdata = simdata, rlatent = err)
}
