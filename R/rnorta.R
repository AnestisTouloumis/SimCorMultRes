rnorta <- function(R = R, cor.matrix = cor.matrix, distr = distr, qparameters = NULL) {
    if (all.equal(R, as.integer(R)) != TRUE | R < 1) 
        stop("'R' must be a positive integer")
    QuantileFunctions <- as.character(distr)
    ans <- rsmvnorm(R = R, cor.matrix = cor.matrix)
    if (length(QuantileFunctions) != ncol(cor.matrix)) 
        stop("'distr' must be a ", ncol(cor.matrix), "-variate vector of strings naming a valid quantile function")
    if (!is.null(qparameters)) {
        qparameters <- as.list(qparameters)
        if (length(qparameters) != ncol(cor.matrix)) 
            stop("'qparameters' must be provided as a list of length ", 
                ncol(cor.matrix))
    }
    ans <- pnorm(ans)
    for (i in seq_len(ncol(cor.matrix))) {
        QuantileFunction <- get(QuantileFunctions[i], mode = "function")
        if (!is.function(QuantileFunction)) 
            stop("Character string ", i, " in `distr' does not correspond to a valid function")
        if (!is.null(qparameters)) 
            formals(QuantileFunction)[pmatch(names(qparameters[[i]]), formalArgs(QuantileFunction))] <- qparameters[[i]]
        ans[, i] <- QuantileFunction(ans[, i])
    }
    ans
}
