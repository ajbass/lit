
# Cross product function
# test_that("Test cross product cpp function", {
#   y1 <- rnorm(10)
#   y2 <- rnorm(10)
#   expect_equal(as.matrix(y1*y2), lit:::pairwise_prod(x = cbind(y1,y2)))
# })

# SVD function
# test_that("Test svd cpp function", {
#   y1 <- rnorm(10)
#   y2 <- rnorm(10)
#   X <- matrix(y1*y2, ncol=2)
#   lit_svd <- lit:::quick_svd(X = X)
#   svd_out <- svd(t(scale(X)))
#   expect_equal(svd_out$d^2, as.numeric(lit_svd$d))
#   expect_equal(sum(abs(lit_svd$U)), sum(abs(svd_out$v)))
# })

# functions in linreg.cpp
test_that("Test quick_lm cpp function", {
  X <- matrix(rnorm(50*5), ncol = 5)
  y <- matrix(rnorm(50*3), ncol = 3)

  lmout <- lit:::.quick_lm_cpp(X = X, Y=y)
  Rlmout <- residuals(lm(y~-1+X))
  expect_equal(as.numeric(Rlmout), as.numeric(lmout))
})

# test_that("Test quick_geno cpp function", {
#   X <- rbinom(n = 100, size =2, prob = 0.25)
#   y <- cbind(1,matrix(rnorm(100*4), ncol = 4))
#   lmout <- lit:::quick_geno(X = y, Y=X)
#   Rlmout <- fitted(lm(X/2~-1 + y))
#   low = 1 / 200
#   high = 1 - low
#   Rlmout[Rlmout<low] <- low
#   Rlmout[Rlmout>high] <- high
#   res = (X - 2 * Rlmout) / sqrt(2 * Rlmout * (1-Rlmout))
#   expect_equal(as.numeric(res), as.numeric(lmout))
# })

# test_that("Test marginal_test function", {
#   X <- matrix(rnorm(50*1), ncol = 1)
#   y <- matrix(rnorm(50*1), ncol = 1)
#   p <- lit:::.marginal_test_F(X = X, y=y, H = matrix(1, nrow(y)))
#   Rlmout <- (lm(y~1 + X))
#   plm <- anova(Rlmout)$`Pr(>F)`[1]#summary(Rlmout)$coef[4]
#   expect_equal(as.numeric(p), as.numeric(plm))
# })

# functions in hsic_cpp
# test_that("Test hsic function", {
#   X <- matrix(rnorm(100 * 100), ncol = 100, nrow = 100)
#   Y <- matrix(rnorm(100 * 100), ncol = 100, nrow = 100)
#   svdX <- svd(X)
#   svdY <- svd(Y)
#   p <- lit:::.hsic_cpp(Kx_v = svdX$v,
#                       Kx_e = svdX$d ^ 2,
#                       Ky_v = svdY$v,
#                       Ky_e = svdY$d ^ 2)
#
#   Xc <- t(X) %*% (X)
#   Yc <- t(Y) %*% (Y)
#   library(CompQuadForm)
#   m = nrow(Y) # number of subjects in study
#   GAMuT = (1/m) * sum(sum(t(Yc) * Xc))
#   Z <- tcrossprod(as.matrix(svdY$d^2), as.matrix(svdX$d^2))
#   Zsort <- sort(Z, decreasing=T)
#   ## derive p-value of GAMuT statistic:
#   scoredavies = GAMuT * m^2
#   results_score <- CompQuadForm::davies(scoredavies, Zsort,acc = 0.000001)
#   davies_pvalue <- (results_score$Qq)
#
#   expect_equal(as.numeric(p), as.numeric(davies_pvalue))
# })

# functions in read_bed.cpp

# functions in lit_cpp.cpp

# check inputs in lit function, and lit_plink
test_that("Test lit function", {
  library(genio)
 # X <- matrix(rbinom(n = 10*10, size =2, prob = 0.25), ncol = 10)
 # y <- matrix(rnorm(10*4), ncol = 4)
 # out_p <- lit::lit(y = y, x = X)

  file <- system.file( "extdata", 'sample.bed', package = "genio", mustWork = TRUE )
  h <- as.matrix(rnorm(10))
  adjustment = as.matrix(rnorm(10))
  # read genotypes and annotation tables
  plink_data <- read_plink( file )
  # genotypes
  X <- plink_data$X
  y <- matrix(rnorm(10*4), ncol = 4)
  out_p <- lit::lit(y = y, x = t(X), adjustment = adjustment, pop_struct = h)

  file <- system.file( "extdata", 'sample.bed', package = "genio", mustWork = TRUE )
  litbed <- lit::lit_plink(y = y, file = file,  adjustment = adjustment, pop_struct=h)
  # Load genio file
  expect_equal(litbed[,7:9], out_p)
})
