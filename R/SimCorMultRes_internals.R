check_cluster_size <- function(cluster_size) {
    if (all.equal(cluster_size, as.integer(cluster_size)) != TRUE |
        cluster_size < 2)
        stop("'clsize' must be a positive integer greater than or equal to two")
}

check_ncategories <- function(categories_no) {
    if (!is.numeric(categories_no) | categories_no < 3)
        stop("'ncategories' must be numeric greater or equal to three")
    if (all.equal(categories_no, as.integer(categories_no)) != TRUE)
        stop("'ncategories' must be a positive integer")
    categories_no
}

check_correlation_matrix <- function(correlation_matrix, cluster_size, rfctn,
    categories_no = NULL) {
    if (!is.numeric(correlation_matrix))
        stop("'cor_matrix' must be numeric")
    if (!is.matrix(correlation_matrix))
        stop("'cor_matrix' must be a matrix")
    if (rfctn == "rbin" | rfctn == "rmult.clm") {
        if (ncol(correlation_matrix) != cluster_size |
            nrow(correlation_matrix) != cluster_size)
            stop("'cor.matrix' must be a ", cluster_size, "x",
                 cluster_size, " matrix")
    } else {
        ncateg2 <- ifelse(rfctn == "rmult.bcl", categories_no,
                          categories_no - 1)
        dimcor <- cluster_size * ncateg2
        if (ncol(correlation_matrix) != dimcor |
            nrow(correlation_matrix) != dimcor)
            stop("'cor_matrix' must be a ", dimcor, "x", dimcor, " matrix")
        for (i in 1:cluster_size) {
            diag_index <- 1:ncateg2 + (i - 1) * ncateg2
            correlation_matrix[diag_index, diag_index] <- diag(1, ncateg2)
        }
    }
    if (!isSymmetric(correlation_matrix))
        stop("'cor_matrix' must be a symmetric matrix")
    if (any(diag(correlation_matrix) != 1))
        stop("the diagonal elements of 'cor_matrix' must be one")
    if (any(correlation_matrix > 1) | any(correlation_matrix < -1))
        stop("all the elements of 'cor_matrix' must be on [-1,1]")
    correlation_matrix_eigen <- eigen(correlation_matrix, symmetric = TRUE,
                                      only.values = TRUE)
    if (any(correlation_matrix_eigen$values <= 0))
        stop("'cor_matrix' must be positive definite")
    correlation_matrix
}

check_xformula <- function(xformula) {
    linear_predictor_formula <- as.formula(xformula)
    if (length(paste0(attr(terms(linear_predictor_formula), "variables"))) == 1)
        stop("No covariates were found in 'formula' ")
    if (attr(terms(linear_predictor_formula), "intercept") == 0) {
        linear_predictor_formula <- update(linear_predictor_formula, ~. + 1)
    }
    linear_predictor_formula
}

check_intercepts <- function(intercepts, cluster_size, rfctn,
                             sample_size = NULL) {
    if (!is.numeric(intercepts))
        stop("'intercepts' must be numeric")
    if (any(is.infinite(intercepts)))
        stop("'intercepts' must not be Inf or -Inf")
    if (rfctn == "rbin") {
        if (!(is.vector(intercepts) & !is.list(intercepts)))
            stop("'intercepts' must be a vector")
        if (length(intercepts) == 1)
            intercepts <- rep(intercepts, cluster_size)
        if (length(intercepts) != cluster_size)
            stop("'intercepts' must have either one or ", cluster_size,
                 " elements")
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
            categories_no <- length(intercepts) + 1
            if (rfctn == "rmult.clm") {
                intercepts <- matrix(intercepts, cluster_size, categories_no -
                  1, TRUE)
                intercepts <- cbind(-Inf, intercepts, Inf)
            } else if (rfctn == "rmult.crm") {
                intercepts <- matrix(intercepts, sample_size,
                                     cluster_size * (categories_no - 1), TRUE)
            } else {
                intercepts <- matrix(intercepts, cluster_size,
                                     categories_no - 1, byrow = TRUE)
            }
        } else {
            categories_no <- ncol(intercepts) + 1
            if (rfctn == "rmult.clm") {
                intercepts <- cbind(-Inf, intercepts, Inf)
                for (i in 1:cluster_size) {
                  if (any(diff(intercepts[i, ]) <= 0))
                    stop("'intercepts' must be increasing at each row")
                }
            } else if (rfctn == "rmult.crm") {
                categories_no <- ncol(intercepts) + 1
                intercepts <- matrix(t(intercepts), sample_size,
                                     cluster_size * (categories_no - 1), TRUE)
            } else {
                categories_no <- ncol(intercepts) + 1
                intercepts <- matrix(intercepts, cluster_size,
                                     categories_no - 1, TRUE)
            }
        }
    }
    intercepts
}

check_betas <- function(betas, cluster_size) {
    if (!(is.vector(betas) & !is.list(betas)) & !is.matrix(betas))
        stop("'betas' must be a vector or a matrix")
    if (!is.numeric(betas))
        stop("'betas' must be numeric")
    if (is.vector(betas) & !is.list(betas)) {
        betas <- rep(betas, cluster_size)
    } else {
        if (nrow(betas) != cluster_size)
            stop("The number of rows in 'betas' should be equal to 'clsize'")
        betas <- c(t(betas))
    }
    betas
}

create_linear_predictor <- function(betas, cluster_size,
                                    linear_predictor_formula, xdata, rfctn,
                                    categories_no = NULL) {
    xmat <- model.matrix(linear_predictor_formula, data = xdata)
    if (rfctn == "rmult.bcl") {
        xmat <- apply(xmat, 2, function(x) rep(x, each = categories_no))
        if (length(betas) != (cluster_size * categories_no * ncol(xmat)))
            stop("The length of 'betas' does not match with the provided covariates") # nolint
    } else {
        xmat <- matrix(xmat[, -1], ncol = ncol(xmat) - 1)
        if (length(betas) != (cluster_size) * ncol(xmat))
            stop("The length of 'betas' does not match with the number of covariates") # nolint
    }
    linear_predictor <- matrix(betas, nrow = nrow(xmat), ncol = ncol(xmat),
        byrow = TRUE) * xmat
    if (rfctn == "rmult.bcl") {
        linear_predictor <- matrix(rowSums(linear_predictor),
                                   ncol = categories_no * cluster_size,
                                   byrow = TRUE)
    } else {
        linear_predictor <- matrix(rowSums(linear_predictor),
                                   ncol = cluster_size, byrow = TRUE)
    }
    as.matrix(linear_predictor)
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

create_rlatent <- function(simulated_latent_variables, sample_size, link,
                           cluster_size, correlation_matrix,
                           rfctn, categories_no = NULL) {
    if (is.null(simulated_latent_variables)) {
        distr <- create_distribution(link)
        correlation_matrix <- check_correlation_matrix(correlation_matrix,
                                                       cluster_size, rfctn,
                                                       categories_no)
        simulated_latent_variables <- rnorta(sample_size, correlation_matrix,
                                             rep(distr,
                                                 nrow(correlation_matrix)))
        if (distr == "qgumbel" & rfctn != "rmult.bcl")
            simulated_latent_variables <- -simulated_latent_variables
    } else {
        if (!is.matrix(simulated_latent_variables))
            stop("'rlatent' must be a matrix")
        if (!is.numeric(simulated_latent_variables))
            stop("'rlatent' must be numeric")
        if (rfctn == "rbin" | rfctn == "rmult.clm") {
            ncol_rlatent <- cluster_size
        } else if (rfctn == "rmult.bcl") {
            ncol_rlatent <- cluster_size * categories_no
        } else {
            ncol_rlatent <- cluster_size * (categories_no - 1)
        }
        if (nrow(simulated_latent_variables) != sample_size |
            ncol(simulated_latent_variables) != ncol_rlatent)
            stop("'rlatent' must be a ", sample_size, "x", ncol_rlatent, " matrix") # nolint
        correlation_matrix <- NULL
        simulated_latent_variables <- simulated_latent_variables
    }
    simulated_latent_variables
}

create_output <- function(simulated_responses, sample_size, cluster_size,
                          simulated_latent_variables, linear_predictor_formula,
                          xdata, rfctn, categories_no = NULL) {
    y <- c(t(simulated_responses))
    id <- rep(1:sample_size, each = cluster_size)
    time <- rep(1:cluster_size, sample_size)
    rownames(simulated_responses) <- rownames(simulated_latent_variables) <-
        paste("i", 1:sample_size, sep = "=")
    colnames(simulated_responses) <- paste("t", 1:cluster_size, sep = "=")
    colnames(simulated_latent_variables) <-
        if (rfctn == "rbin" | rfctn == "rmult.clm") {
            paste("t", 1:cluster_size, sep = "=")
            } else if (rfctn == "rmult.bcl") {
                paste("t=", rep(1:cluster_size, each = categories_no), " & j=",
                      rep(1:categories_no, cluster_size), sep = "")
                } else {
                    paste("t=", rep(1:cluster_size, each = categories_no - 1),
                          " & j=", rep(1:(categories_no - 1), cluster_size),
                          sep = "")
    }
    sim_model_frame <- model.frame(formula = linear_predictor_formula,
                                   data = xdata)
    simdata <- data.frame(y, sim_model_frame, id, time)
    list(Ysim = simulated_responses, simdata = simdata,
         rlatent = simulated_latent_variables)
}


apply_threshold <- function(linear_predictor, simulated_latent_variables,
                            cluster_size, rfctn, intercepts = NULL,
                            categories_no = NULL) {
    sample_size <- nrow(linear_predictor)
    if (rfctn == "rmult.clm" | rfctn == "rmult.crm") {
        u_sim <- simulated_latent_variables - linear_predictor
    } else {
        u_sim <- simulated_latent_variables + linear_predictor
    }
    if (rfctn == "rbin") {
        simulated_responses <- matrix(0, sample_size, cluster_size)
        for (i in 1:cluster_size) simulated_responses[, i] <-
                cut(u_sim[, i] - 2 * linear_predictor[, i], intercepts[i, ],
                    labels = FALSE)
        simulated_responses <- 2 - simulated_responses
    } else if (rfctn == "rmult.bcl") {
        u_sim <- matrix(as.vector(t(u_sim)), nrow = cluster_size * sample_size,
                        ncol = categories_no, byrow = TRUE)
        simulated_responses <- apply(u_sim, 1, which.max)
        simulated_responses <- matrix(simulated_responses, ncol = cluster_size,
                                      byrow = TRUE)
    } else if (rfctn == "rmult.clm") {
        simulated_responses <- matrix(0, sample_size, cluster_size)
        for (i in 1:cluster_size) simulated_responses[, i] <-
                cut(u_sim[, i], intercepts[i, ], labels = FALSE)
    } else {
        simulated_responses <- matrix(as.numeric(t(u_sim <= intercepts)),
                        sample_size * cluster_size, categories_no - 1, TRUE)
        for (i in 1:(categories_no - 1))
            simulated_responses[, i] <- ifelse(simulated_responses[, i] == 1,
                                               i, categories_no)
        simulated_responses <- apply(simulated_responses, 1, min)
        simulated_responses <- matrix(simulated_responses, sample_size,
                                      cluster_size, TRUE)
    }
    simulated_responses
}

create_betas_acl2bcl <- function(intercepts = intercepts,
                                 categories_no = categories_no, betas = betas) {
    intercepts_bcl <- t(apply(intercepts, 1, function(z) rev(cumsum(rev(z)))))
    cluster_size <- nrow(intercepts)
    dim_betas <- length(betas) / cluster_size
    betas_matrix <- matrix(betas, cluster_size, dim_betas, TRUE)
    betas_matrix_bcl <- t(apply(betas_matrix, 1, function(z)
        rep(categories_no - seq(categories_no - 1), each = length(z)) * z))
    betas_bcl <- matrix(0, cluster_size,
                        categories_no - 1 + ncol(betas_matrix_bcl))
    for (i in seq(categories_no - 1)) {
        betas_bcl[, ((i - 1) * (dim_betas + 1) + 1):(i * (dim_betas + 1))] <-
          cbind(intercepts_bcl[, i],
            betas_matrix_bcl[, ((i - 1) * (dim_betas) + 1):(i * (dim_betas))])
    }
    betas_bcl <- cbind(betas_bcl, matrix(0, cluster_size, dim_betas + 1))
    betas_bcl
}
