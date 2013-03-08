x <- rnorm(2^24)
y <- x
X <- cbind(1, x)
system.time(for (i in 1:5) lm.fit(X, y))
