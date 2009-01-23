
twoGaussiansNull <- function(x, p.adj.method="BY", max.adj.p = 0.1, var.equal=FALSE, ...){
    require("mclust")
    ## check arguments:
    stopifnot(is.numeric(x), p.adj.method %in% p.adjust.methods)
    ## set model
    thisModel <- ifelse(var.equal, "E", "V")
    ## fit two gaussians
    xmclust <- Mclust(na.omit(x), G = 2, modelNames = thisModel, ...)
    ## convert probe intensities to p-values under the assumption that
    ##  the lower model is the null model
    nu <- which.min(xmclust$parameters$mean)
    xp <- pnorm(x, xmclust$parameters$mean[nu], sqrt(xmclust$parameters$variance$sigmasq[nu]), lower.tail = FALSE)
    xq <- p.adjust(xp, method = p.adj.method)
    min(x[xq <= max.adj.p], na.rm = TRUE)
}# twoGaussiansNull
