# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

profileloglik <- function(sigma, S, n) {
    .Call('_gicf_profileloglik', PACKAGE = 'gicf', sigma, S, n)
}

gicf_wrapper <- function(start, adj, n, S, lambda, lambda_max, tolout, tolin, iterout, iterin) {
    .Call('_gicf_gicf_wrapper', PACKAGE = 'gicf', start, adj, n, S, lambda, lambda_max, tolout, tolin, iterout, iterin)
}

