library(boot)
for (i in 1 : 2) {
    ## example from boot:

    nuke <- nuclear[, c(1, 2, 5, 7, 8, 10, 11)]
    nuke.lm <- glm(log(cost) ~ date+log(cap)+ne+ct+log(cum.n)+pt, data = nuke)
    nuke.diag <- glm.diag(nuke.lm)
    nuke.res <- nuke.diag$res * nuke.diag$sd
    nuke.res <- nuke.res - mean(nuke.res)
    
    nuke.data <- data.frame(nuke, resid = nuke.res, fit = fitted(nuke.lm))
    
    new.data <- data.frame(cost = 1, date = 73.00, cap = 886, ne = 0,
                           ct = 0, cum.n = 11, pt = 1)
    new.fit <- predict(nuke.lm, new.data)
    
    nuke.fun <- function(dat, inds, i.pred, fit.pred, x.pred)
    {
        lm.b <- glm(fit+resid[inds] ~ date+log(cap)+ne+ct+log(cum.n)+pt,
                    data = dat)
        pred.b <- predict(lm.b, x.pred)
        c(coef(lm.b), pred.b - (fit.pred + dat$resid[i.pred]))
    }
    
    nuke.boot <- boot(nuke.data, nuke.fun, R = 999, m = 1, 
                      fit.pred = new.fit, x.pred = new.data)
    mean(nuke.boot$t[, 8]^2)
    new.fit - sort(nuke.boot$t[, 8])[c(975, 25)]
}
