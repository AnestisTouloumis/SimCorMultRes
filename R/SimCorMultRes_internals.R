check_cluster_size <- function(clsize) {
    if (all.equal(clsize, as.integer(clsize)) != TRUE | clsize < 2)
        stop("'clsize' must be a positive integer greater than or equal to two")
}

check_ncategories <- function(ncategories) {
    if (!is.numeric(ncategories) | ncategories < 3)
        stop("'ncategories' must be numeric greater or equal to three")
    if (all.equal(ncategories, as.integer(ncategories)) != TRUE)
        stop("'ncategories' must be a positive integer")
    ncategories
}

check_correlation_matrix <- function(cor_matrix, clsize, rfctn,
    ncategories = NULL) {
    if (!is.numeric(cor_matrix))
        stop("'cor_matrix' must be numeric")
    if (!is.matrix(cor_matrix))
        stop("'cor_matrix' must be a matrix")
    if (rfctn == "rbin" | rfctn == "rmult.clm") {
        if (ncol(cor_matrix) != clsize | nrow(cor_matrix) != clsize)
            stop("'cor_matrix' must be a ", clsize, "x", clsize, " matrix")
    } else {
        ncateg2 <- ifelse(rfctn == "rmult.bcl", ncategories, ncategories - 1)
        dimcor <- clsize * ncateg2
        if (ncol(cor_matrix) != dimcor | nrow(cor_matrix) != dimcor)
            stop("'cor_matrix' must be a ", dimcor, "x", dimcor, " matrix")
        for (i in 1:clsize) {
            diag_index <- 1:ncateg2 + (i - 1) * ncateg2
            cor_matrix[diag_index, diag_index] <- diag(1, ncateg2)
        }
    }
    if (!isSymmetric(cor_matrix))
        stop("'cor_matrix' must be a symmetric matrix")
    if (any(diag(cor_matrix) != 1))
        stop("the diagonal elements of 'cor_matrix' must be one")
    if (any(cor_matrix > 1) | any(cor_matrix < -1))
        stop("all the elements of 'cor_matrix' must be on [-1,1]")
    if (any(eigen(cor_matrix, symmetric = TRUE, only.values = TRUE)$values <=
        0))
        stop("'cor_matrix' must be positive definite")
    cor_matrix
}

check_xformula <- function(xformula) {
    lpformula <- as.formula(xformula)
    if (length(paste0(attr(terms(lpformula), "variables"))) == 1)
        stop("No covariates were found in 'formula' ")
    if (attr(terms(lpformula), "intercept") == 0) {
        lpformula <- update(lpformula, ~. + 1)
    }
    lpformula
}

check_intercepts <- function(intercepts, clsize, rfctn, sample_size = NULL) {
    if (!is.numeric(intercepts))
        stop("'intercepts' must be numeric")
    if (any(is.infinite(intercepts)))
        stop("'intercepts' must not be Inf or -Inf")
    if (rfctn == "rbin") {
        if (!(is.vector(intercepts) & !is.list(intercepts)))
            stop("'intercepts' must be a vector")
        if (length(intercepts) == 1)
            intercepts <- rep(intercepts, clsize)
        if (length(intercepts) != clsize)
            stop("'intercepts' must have either one or ", clsize, " elements")
        intercepts <- cbind(-Inf, intercepts, Inf)
    } else {
        if (!(is.vector(intercepts) & !is.list(intercepts)) &
            !is.matrix(intercepts))
            stop("'intercepts' must be a vector or a matrix")
        if (is.vector(intercepts) & !is.list(intercepts)) {
            if (length(intercepts) == 1)
                stop("'intercepts' must have at least 2 elements")
            if (rfctn == "rmult.clm" & any(diff(intercepts) <= 0))
                stop("'intercepts' must be increasing")
            ncategories <- length(intercepts) + 1
            if (rfctn == "rmult.clm") {
                intercepts <- matrix(intercepts, clsize, ncategories -
                  1, TRUE)
                intercepts <- cbind(-Inf, intercepts, Inf)
            } else if (rfctn == "rmult.crm") {
                intercepts <- matrix(intercepts, sample_size,
                                     clsize * (ncategories - 1), byrow = TRUE)
            } else {
                intercepts <- matrix(intercepts, clsize, ncategories -
                  1, byrow = TRUE)
            }
        } else {
            ncategories <- ncol(intercepts) + 1
            if (rfctn == "rmult.clm") {
                intercepts <- cbind(-Inf, intercepts, Inf)
                for (i in 1:clsize) {
                  if (any(diff(intercepts[i, ]) <= 0))
                    stop("'intercepts' must be increasing at each row")
                }
            } else if (rfctn == "rmult.crm") {
                ncategories <- ncol(intercepts) + 1
                intercepts <- matrix(t(intercepts), sample_size,
                                     clsize * (ncategories - 1), byrow = TRUE)
            } else {
                ncategories <- ncol(intercepts) + 1
                intercepts <- matrix(intercepts, clsize, ncategories -
                  1, byrow = TRUE)
            }
        }
    }
    intercepts
}

check_betas <- function(betas, clsize) {
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

create_linear_predictor <- function(betas, clsize, lpformula, xdata, rfctn,
    ncategories = NULL) {
    xmat <- model.matrix(lpformula, data = xdata)
    if (rfctn == "rmult.bcl") {
        xmat <- apply(xmat, 2, function(x) rep(x, each = ncategories))
        if (length(betas) != (clsize * ncategories * ncol(xmat)))
            stop("The length of 'betas' does not match with the provided covariates") # nolint
    } else {
        xmat <- matrix(xmat[, -1], ncol = ncol(xmat) - 1)
        if (length(betas) != (clsize) * ncol(xmat))
            stop("The length of 'betas' does not match with the number of covariates") # nolint
    }
    lin_pred <- matrix(betas, nrow = nrow(xmat), ncol = ncol(xmat),
        byrow = TRUE) * xmat
    if (rfctn == "rmult.bcl") {
        lin_pred <- matrix(rowSums(lin_pred), ncol = ncategories * clsize,
            byrow = TRUE)
    } else {
        lin_pred <- matrix(rowSums(lin_pred), ncol = clsize, byrow = TRUE)
    }
    as.matrix(lin_pred)
}

create_distribution <- function(link) {
    if (length(link) != 1)
        stop("The length of 'link' must be one")
    links <- c("probit", "logit", "cloglog", "cauchit")
    if (!is.element(link, links))
        stop("'link' must be 'probit','logit','cloglog' and/or 'cauchit'")
    distr <- switch(link, probit = "qnorm", logit = "qlogis",
        cloglog = "qgumbel", cauchit = "qcauchy")
    distr
}

create_rlatent <- function(rlatent, sample_size, link, clsize, cor_matrix,
                           rfctn, ncategories = NULL) {
    if (is.null(rlatent)) {
        distr <- create_distribution(link)
        cor_matrix <- check_correlation_matrix(cor_matrix, clsize, rfctn,
            ncategories)
        rlatent <- rnorta(sample_size, cor_matrix, rep(distr, nrow(cor_matrix)))
        if (distr == "qgumbel" & rfctn != "rmult.bcl")
            rlatent <- -rlatent
    } else {
        if (!is.matrix(rlatent))
            stop("'rlatent' must be a matrix")
        if (!is.numeric(rlatent))
            stop("'rlatent' must be numeric")
        if (rfctn == "rbin" | rfctn == "rmult.clm") {
            ncol_rlatent <- clsize
        } else if (rfctn == "rmult.bcl") {
            ncol_rlatent <- clsize * ncategories
        } else {
            ncol_rlatent <- clsize * (ncategories - 1)
        }
        if (nrow(rlatent) != sample_size | ncol(rlatent) != ncol_rlatent)
            stop("'rlatent' must be a ", sample_size, "x", ncol_rlatent, " matrix") # nolint
        cor_matrix <- NULL
        rlatent <- rlatent
    }
    rlatent
}

create_output <- function(y_sim, sample_size, clsize, rlatent, lpformula, xdata,
                          rfctn, ncategories = NULL) {
    y <- c(t(y_sim))
    id <- rep(1:sample_size, each = clsize)
    time <- rep(1:clsize, sample_size)
    rownames(y_sim) <- rownames(rlatent) <- paste("i", 1:sample_size, sep = "=")
    colnames(y_sim) <- paste("t", 1:clsize, sep = "=")
    colnames(rlatent) <- if (rfctn == "rbin" | rfctn == "rmult.clm") {
        paste("t", 1:clsize, sep = "=")
    } else if (rfctn == "rmult.bcl") {
        paste("t=", rep(1:clsize, each = ncategories), " & j=",
              rep(1:ncategories, clsize), sep = "")
    } else {
        paste("t=", rep(1:clsize, each = ncategories - 1), " & j=",
            rep(1:(ncategories - 1), clsize), sep = "")
    }
    sim_model_frame <- model.frame(formula = lpformula, data = xdata)
    simdata <- data.frame(y, sim_model_frame, id, time)
    list(Ysim = y_sim, simdata = simdata, rlatent = rlatent)
}


apply_threshold <- function(lin_pred, rlatent, clsize, rfctn, intercepts = NULL,
    ncategories = NULL) {
    sample_size <- nrow(lin_pred)
    if (rfctn == "rmult.clm" | rfctn == "rmult.crm") {
        u_sim <- rlatent - lin_pred
    } else {
        u_sim <- rlatent + lin_pred
    }
    if (rfctn == "rbin") {
        y_sim <- matrix(0, sample_size, clsize)
        for (i in 1:clsize) y_sim[, i] <- cut(u_sim[, i] - 2 * lin_pred[, i],
            intercepts[i, ], labels = FALSE)
        y_sim <- 2 - y_sim
    } else if (rfctn == "rmult.bcl") {
        u_sim <- matrix(as.vector(t(u_sim)), nrow = clsize * sample_size,
                        ncol = ncategories, byrow = TRUE)
        y_sim <- apply(u_sim, 1, which.max)
        y_sim <- matrix(y_sim, ncol = clsize, byrow = TRUE)
    } else if (rfctn == "rmult.clm") {
        y_sim <- matrix(0, sample_size, clsize)
        for (i in 1:clsize) y_sim[, i] <- cut(u_sim[, i], intercepts[i, ],
            labels = FALSE)
    } else {
        y_sim <- matrix(as.numeric(t(u_sim <= intercepts)),
                        sample_size * clsize, ncategories - 1, TRUE)
        for (i in 1:(ncategories - 1)) y_sim[, i] <- ifelse(y_sim[, i] ==
            1, i, ncategories)
        y_sim <- apply(y_sim, 1, min)
        y_sim <- matrix(y_sim, sample_size, clsize, byrow = TRUE)
    }
    y_sim
}

create_betas_acl2bcl <- function(intercepts = intercepts,
    ncategories = ncategories, betas = betas) {
    intercepts_bcl <- t(apply(intercepts, 1, function(z) rev(cumsum(rev(z)))))
    clsize <- nrow(intercepts)
    dim_betas <- length(betas) / clsize
    betas_matrix <- matrix(betas, clsize, dim_betas, TRUE)
    betas_matrix_bcl <- t(apply(betas_matrix, 1, function(z) rep(ncategories -
        seq(ncategories - 1), each = length(z)) * z))
    betas_bcl <- matrix(0, clsize, ncategories - 1 + ncol(betas_matrix_bcl))
    for (i in seq(ncategories - 1)) {
        betas_bcl[, ((i - 1) * (dim_betas + 1) + 1):(i * (dim_betas + 1))] <-
          cbind(intercepts_bcl[, i],
            betas_matrix_bcl[, ((i - 1) * (dim_betas) + 1):(i * (dim_betas))])
    }
    betas_bcl <- cbind(betas_bcl, matrix(0, clsize, dim_betas + 1))
    betas_bcl
}
