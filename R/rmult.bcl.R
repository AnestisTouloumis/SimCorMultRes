rmult.bcl <- function(clsize = clsize, ncategories = ncategories, betas = betas, xformula = formula(xdata), xdata = parent.frame(), cor.matrix = cor.matrix, 
    rlatent = NULL) {
    if (all.equal(clsize, as.integer(clsize)) != TRUE | clsize < 2) 
        stop("'clsize' must be a positive integer greater than or equal to two")
    if (!is.numeric(ncategories) | ncategories < 3) 
        stop("'ncategories' must be greater than or equal to three")
    ncategories <- as.integer(ncategories)
    if (all.equal(ncategories, as.integer(ncategories)) != TRUE | ncategories < 3) 
        stop("'ncategories' must be a positive integer greater than or equal to three")
    if (!is.vector(betas) & !is.matrix(betas)) 
        stop("'betas' must be a vector or a matrix")
    if (!is.numeric(betas)) 
        stop("'betas' must be numeric")
    if (is.vector(betas)) {
        betas <- rep(betas, clsize)
    } else {
        if (nrow(betas) != clsize) 
            stop("The number of rows in 'betas' should be equal to 'clsize'")
        betas <- c(t(betas))
    }
    lpformula <- as.formula(xformula)
    if (attr(terms(lpformula), "intercept") == 0) 
        stop("'formula' must have an intercept term")
    Xmat <- model.matrix(lpformula, data = xdata)
    Xmat <- apply(Xmat, 2, function(x) rep(x, each = ncategories))
    if (length(betas) != (clsize * ncategories * ncol(Xmat))) 
        stop("The length of 'betas' does not match with the provided covariates")
    lin.pred <- matrix(betas, nrow = nrow(Xmat), ncol = ncol(Xmat), byrow = TRUE) * Xmat
    lin.pred <- matrix(rowSums(lin.pred), ncol = ncategories * clsize, byrow = TRUE)
    ncol.lp <- clsize * ncategories
    R <- nrow(lin.pred)
    if (is.null(rlatent)) {
        if (!is.matrix(cor.matrix)) 
            stop("'cor.matrix' must be a matrix")
        if (!is.numeric(cor.matrix)) 
            stop("'cor.matrix' must be numeric")
        if (ncol(cor.matrix) != ncol.lp | nrow(cor.matrix) != ncol.lp) 
            stop("'cor.matrix' must be a ", ncol.lp, "x", ncol.lp, " matrix")
        if (!isSymmetric(cor.matrix)) 
            stop("'cor.matrix' must be a symmetric matrix")
        for (i in 1:clsize) {
            diag.index <- 1:ncategories + (i - 1) * ncategories
            cor.matrix[diag.index, diag.index] <- diag(1, ncategories)
        }
        if (any(cor.matrix > 1) | any(cor.matrix < -1)) 
            stop("all the elements of 'cor.matrix' must be on [-1,1]")
        if (any(eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <= 0)) 
            stop("'cor.matrix' must respect the local independence of the alternatives and must be positive definite")
        err <- rnorta(R, cor.matrix, rep("qgumbel", ncol.lp))
    } else {
        if (!is.matrix(rlatent)) 
            stop("'rlatent' must be a matrix")
        if (!is.numeric(rlatent)) 
            stop("'rlatent' must be numeric")
        if (ncol(rlatent) != ncol.lp | nrow(rlatent) != R) 
            stop("'rlatent' must be a ", R, "x", ncol.lp, " matrix")
        cor.matrix <- NULL
        err <- rlatent
    }
    U <- lin.pred + err
    U <- matrix(as.vector(t(U)), nrow = clsize * R, ncol = ncategories, TRUE)
    Ysim <- apply(U, 1, which.max)
    Ysim <- matrix(Ysim, ncol = clsize, byrow = TRUE)
    id <- rep(1:R, each = clsize)
    time <- rep(1:clsize, R)
    y <- c(t(Ysim))
    lpformula <- update(lpformula, ~. - 1)
    rownames(Ysim) <- rownames(err) <- paste("i", 1:R, sep = "=")
    colnames(Ysim) <- paste("t", 1:clsize, sep = "=")
    colnames(err) <- paste("t=", rep(1:clsize, each = ncategories), " & j=", rep(1:ncategories, clsize), sep = "")
    simdata <- data.frame(y, model.frame(formula = lpformula, data = xdata), id, time)
    list(Ysim = Ysim, simdata = simdata, rlatent = err)
}
