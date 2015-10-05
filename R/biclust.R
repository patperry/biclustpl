#  Copyright 2015 Cheryl J. Flynn and Patrick O. Perry.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.


validate_clusters <- function(n, clusters)
{
    if (anyNA(clusters))
        stop("missing value in 'clusters' argument")

    if (is.factor(clusters)) {
        nclusters <- nlevels(clusters)
        clusters <- as.integer(clusters)
    } else {
        clusters <- as.integer(clusters)
        if (length(clusters) == 1) {
            nclusters <- clusters
            clusters <- sample.int(nclusters, n, replace=TRUE)
        } else {
            nclusters <- max(clusters)
        }

        if (nclusters < 1)
            stop("number of clusters must be >= 1")
    }

    if (length(clusters) != n) {
        stop("'clusters' argument must have length '", n, "'")
    }

    list(nclusters=nclusters, clusters=clusters)
}


biclust_dense <- function(x, row_clusters, col_clusters,
                          family = c("gaussian", "poisson", "binomial"),
                          nstart = 1,
                          epsilon = 1e-6, maxit = 10000L, trace = FALSE)
{
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    m <- nrow(x)
    n <- ncol(x)

    row <- validate_clusters(nrow(x), row_clusters)
    col <- validate_clusters(ncol(x), col_clusters)

    family <- match.arg(family)

    nstart <- round(nstart)
    if (nstart < 1)
        stop("invalid 'nstart' argument")

    storage.mode(epsilon) <- "double"
    storage.mode(maxit) <- "integer"

    K <- row$nclusters
    L <- col$nclusters
    row_cl0 <- row$clusters - 1L
    col_cl0 <- col$clusters - 1L

    best <- NULL

    for (start in seq_len(nstart)) {
        ans <- .Call(C_biclust_dense, x, K, row_cl0, L, col_cl0,
                     family, epsilon, maxit, trace)
        if (start == 1 || ans$loglik > best$loglik) {
            best <- ans
        }
        if (start < nstart) {
            row_cl0 <- sample.int(K, m, replace=TRUE) - 1L
            col_cl0 <- sample.int(L, n, replace=TRUE) - 1L
        }
    }

    ans <- best

    ans$row_clusters <- ans$row_clusters + 1L
    ans$col_clusters <- ans$col_clusters + 1L

    ans$means <- ans$sums / ans$sizes
    ans$means[ans$sizes == 0] <- NA

    ans$row_means <- ans$row_sums / ans$row_sizes
    ans$row_means[ans$row_sizes == 0] <- NA

    ans$col_means <- ans$col_sums / ans$col_sizes
    ans$col_means[ans$col_sizes == 0] <- NA

    ans
}
