rmult.crm <- function(clsize = clsize, intercepts = intercepts, betas = betas, 
    xformula = formula(xdata), xdata = parent.frame(), link = "logit", 
    cor.matrix = cor.matrix, rlatent = NULL) {
    if (all.equal(clsize, as.integer(clsize)) != TRUE | clsize < 2) 
        stop("'clsize' must be a positive integer greater than or equal to two")
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
    if (length(paste0(attr(terms(lpformula), "variables"))) == 1) 
        stop("No covariates were found in 'formula' ")
    Xmat <- model.matrix(lpformula, data = xdata)
    if (attr(terms(lpformula), "intercept") == 1) {
        Xmat <- matrix(Xmat[, -1], ncol = ncol(Xmat) - 1)
    } else {
        att <- attr(Xmat, "assign")
        factor.columns <- unique(att[duplicated(att)])
        if (length(factor.columns) >= 1) {
            Xmat <- Xmat[, -match(factor.columns[1], att)]
            Xmat <- data.matrix(Xmat)
        }
    }
    if (length(betas) != (clsize) * ncol(Xmat)) 
        stop("The length of 'betas' does not match with the number of covariates")
    lin.pred <- matrix(betas, nrow = nrow(Xmat), ncol = ncol(Xmat), byrow = TRUE) * 
        Xmat
    lin.pred <- matrix(rowSums(lin.pred), ncol = clsize, byrow = TRUE)
    lin.pred <- as.matrix(lin.pred)
    R <- nrow(lin.pred)
    if (!is.vector(intercepts) & !is.matrix(intercepts)) 
        stop("'intercepts' must be a vector or a matrix")
    if (!is.numeric(intercepts)) 
        stop("'intercepts' must be numeric")
    if (is.vector(intercepts)) {
        if (length(intercepts) == 1) 
            stop("'intercepts' must have at least 2 elements")
        ncategories <- length(intercepts) + 1
        intercepts <- matrix(intercepts, R, clsize * (ncategories - 1), 
            byrow = TRUE)
    } else {
        ncategories <- ncol(intercepts) + 1
        intercepts <- matrix(t(intercepts), R, clsize * (ncategories - 
            1), byrow = TRUE)
    }
    lin.pred.extended <- t(apply(lin.pred, 1, function(x) rep(x, each = ncategories - 
        1)))
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
        if (!is.numeric(cor.matrix)) 
            stop("'cor.matrix' must be numeric")
        cor.matrix <- as.matrix(cor.matrix)
        dimcor <- clsize * (ncategories - 1)
        if (ncol(cor.matrix) != dimcor) 
            stop("'cor.matrix' must be a ", dimcor, "x", dimcor, " matrix")
        if (!isSymmetric(cor.matrix)) 
            stop("'cor.matrix' must be symmetric")
        for (i in 1:clsize) {
            diag.index <- 1:(ncategories - 1) + (i - 1) * (ncategories - 
                1)
            cor.matrix[diag.index, diag.index] <- diag(1, ncategories - 
                1)
        }
        if (any(cor.matrix > 1) | any(cor.matrix < -1)) 
            stop("all the elements of 'cor.matrix' must be on [-1,1]")
        if (any(eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <= 
            0)) 
            stop("'cor.matrix' must be positive definite")
        err <- rnorta(R, cor.matrix, rep(distr, ncol(cor.matrix)))
        index <- distr == "qgumbel"
        if (distr == "qgumbel") 
            err <- -err
    } else {
        if (!is.matrix(rlatent)) 
            stop("'rlatent' must be a matrix")
        if (!is.numeric(rlatent)) 
            stop("'rlatent' must be numeric")
        if (ncol(rlatent) != ncol(lin.pred.extended) | nrow(rlatent) != 
            R) 
            stop("'rlatent' must be a ", R, "x", ncol(lin.pred.extended), 
                " matrix")
        cor.matrix <- NULL
        err <- rlatent
    }
    U <- err - lin.pred.extended
    Ysim <- matrix(as.numeric(t(U <= intercepts)), R * clsize, ncategories - 
        1, TRUE)
    for (i in 1:(ncategories - 1)) Ysim[, i] <- ifelse(Ysim[, i] == 1, 
        i, ncategories)
    Ysim <- apply(Ysim, 1, min)
    Ysim <- matrix(Ysim, R, clsize, byrow = TRUE)
    id <- rep(1:R, each = clsize)
    time <- rep(1:clsize, R)
    y <- c(t(Ysim))
    rownames(Ysim) <- rownames(err) <- paste("i", 1:R, sep = "=")
    colnames(Ysim) <- paste("t", 1:clsize, sep = "=")
    colnames(err) <- paste("t=", rep(1:clsize, each = ncategories - 1), 
        " & j=", rep(1:(ncategories - 1), clsize), sep = "")
    simdata <- data.frame(y, model.frame(formula = lpformula, data = xdata), 
        id, time)
    list(Ysim = Ysim, simdata = simdata, rlatent = err)
}
