#' Latent Interaction Testing
#'
#' @description
#' \code{lit} performs a kernel-based testing procedure, Latent Interaction Testing (LIT), using a set of traits and SNPs.
#' LIT tests whether the squared residuals (SQ) and cross products (CP) are statistically independent of the genotypes.
#' In particular, we construct a kernel matrix for the SQ/CP terms to measure the pairwise
#' similarity between individuals, and also construct an analogous one for the genotypes.
#' We then test whether these two matrices are independent.
#' Currently, we implement the linear and projection kernel functions to measure pairwise similarity between individuals.
#' We then combine the p-values of these implementations using a Cauchy combination test to maximize the number of discoveries.
#'
#' @param y matrix of traits (n observations by k traits)
#' @param x matrix of SNPs (n observations by m SNPs)
#' @param adjustment matrix of covariates to adjust traits
#' @param pop_struct matrix of PCs that captures population structure
#'
#' @return
#' A data frame of p-values where the columns are
#'  \itemize{
#'  \item{\code{wlit}: }{LIT using a linear kernel}
#'  \item{\code{ulit}: }{LIT using a projection kernel}
#'  \item{\code{alit}: }{Cauchy combination test of the above two LIT implementations.}
#' }
#'
#' @examples
#' # set seed
#' set.seed(123)
#'
#' # Generate SNPs and traits
#' X <- matrix(rbinom(10*2, size = 2, prob = 0.25), ncol = 2)
#' Y <- matrix(rnorm(10*4), ncol = 4)
#'
#' out <- lit(Y, X)
#'
#' @seealso \code{\link{lit_plink}}
#' @export
lit <- function(y,
                x,
                adjustment = NULL,
                pop_struct = NULL) {
  if (!(is.matrix(y) & is.matrix(x))) {
    stop("Inputs must be matrices.")
  }
  if (!(nrow(y) == nrow(x))) {
    stop("Observations are rows for each input.")
  }
  if (anyNA(y)) {
    stop("Current implementation requires no NAs in trait matrix. Need to adjust inputs accordingly ... feature will be added soon.")
  }

  if (is.null(pop_struct)) {
    h = matrix(1, nrow = nrow(x), ncol = 1)
  } else {
    if (!(is.matrix(pop_struct) & (nrow(y) == nrow(pop_struct)))){
      stop("pop_struct must be matrix and rows are observations.")
    }
    h = cbind(1, pop_struct)
  }

  if (is.null(adjustment)) {
    y <- .quick_lm_cpp(Xs = h, Ys = y)
  } else {
    if (!(is.matrix(adjustment) & (nrow(y) == nrow(adjustment)))){
      stop("pop_struct must be matrix and rows are observations.")
    }
     y <- .quick_lm_cpp(Xs = cbind(h, adjustment), Ys = y)
  }
  m <- ncol(x)
  p <- matrix(nrow = m, ncol = 3)
  ind <- !is.na(x)
  for (i in 1:m) {
    indtmp <- ind[,i] # remove NAs from genotypes
    xtmp <- x[indtmp, i, drop=F]
    ytmp <- y[indtmp,,drop = F]
    htmp <- h[indtmp,,drop = F]
    p[i,] <- .lit_cpp(Xs = xtmp, Ys = ytmp, Hs = htmp)
  }
  colnames(p) <-  c("wlit", "ulit", "alit")
  as.data.frame(p)
}

#' LIT correcting for dominance effects
#'
#' @description
#' Internal use for now
#'
#' @param y matrix of traits (n observations by k traits)
#' @param x matrix of SNPs (n observations by m SNPs)
#' @param adjustment matrix of covariates to adjust traits
#' @param pop_struct matrix of PCs that captures population structure
lit_h <- function(y,
                  x,
                  adjustment = NULL,
                  pop_struct = NULL) {
  if (!(is.matrix(y) & is.matrix(x))) {
    stop("Inputs must be matrices.")
  }
  if (!(nrow(y) == nrow(x))) {
    stop("Observations are rows for each input.")
  }
  if (anyNA(y)) {
    stop("Current implementation requires no NAs in trait matrix. Need to adjust inputs accordingly ... feature will be added soon.")
  }

  if (is.null(pop_struct)) {
    h = matrix(1, nrow = nrow(x), ncol = 1)
  } else {
    if (!(is.matrix(pop_struct) & (nrow(y) == nrow(pop_struct)))){
      stop("pop_struct must be matrix and rows are observations.")
    }
    h = cbind(1, pop_struct)
  }

  if (is.null(adjustment)) {
      y <- .quick_lm_cpp(Xs = h, Ys = y)
  } else {
    if (!(is.matrix(adjustment) & (nrow(y) == nrow(adjustment)))){
      stop("pop_struct must be matrix and rows are observations.")
    }
    y <- .quick_lm_cpp(Xs = cbind(h, adjustment), Ys = y)
  }
  m <- ncol(x)
  p <- matrix(nrow = m, ncol = 3)
  ind <- !is.na(x)
  for (i in 1:m) {
    indtmp <- ind[,i] # remove NAs from genotypes
    xtmp <- x[indtmp, i, drop=F]
    ytmp <- y[indtmp,,drop = F]
    htmp <- h[indtmp,,drop = F]
    maf <-   mean(xtmp) / 2
    if (maf > 0.5) xtmp = 2 - xtmp
    m1 <- as.numeric(xtmp == 1)
    m2 <- as.numeric(xtmp == 2)
    xdom <- model.matrix(~htmp + m1 + m2)
    ytmp <- residuals(lm(ytmp ~ xdom))
    p[i,] <- .lit_cpp(Xs = xtmp, Ys = ytmp, Hs = htmp)
  }
  colnames(p) <-  c("wlit", "ulit", "alit")
  as.data.frame(p)
}

#' Latent Interaction Testing
#'
#' \code{lit_plink} performs a kernel-based testing procedure, Latent Interaction Testing (LIT), using a set of traits and SNPs.
#' LIT tests whether the squared residuals (SQ) and cross products (CP) are statistically independent of the genotypes.
#' In particular, we construct a kernel matrix for the SQ/CP terms to measure the pairwise
#' similarity between individuals, and also construct an analogous one for the genotypes.
#' We then test whether these two matrices are independent.
#' Currently, we implement the linear and projection kernel functions to measure pairwise similarity between individuals.
#' We then combine the p-values of these implementations using a Cauchy combination test to maximize the number of discoveries.
#' This function is suitable for large  datasets (e.g., UK Biobank) in plink format.
#'  Note that our code to process plink files builds from the
#' \code{\link{genio}} R package
#'
#' @param y matrix of traits (n observations by k traits)
#' @param file path to plink files
#' @param adjustment matrix of covariates to adjust traits
#' @param pop_struct matrix of PCs that captures population structure
#' @param verbose If TRUE (default) print progress.
#'
#' @return
#' A data frame of p-values where the columns are
#'  \itemize{
#'  \item{\code{wlit}: }{LIT using a linear kernel}
#'  \item{\code{ulit}: }{LIT using a projection kernel}
#'  \item{\code{alit}: }{Cauchy combination test of the above two LIT implementations.}
#' }
#'
#' @examples
#' # set seed
#' set.seed(123)
#'
#' # path to plink files
#' file <- system.file("extdata", 'sample.bed', package = "genio", mustWork = TRUE)
#'
#' # Generate trait expression
#' Y <- matrix(rnorm(10*4), ncol = 4)
#'
#' out <- lit_plink(Y, file = file)
#'
#' @seealso \code{\link{lit}}
#' @export
lit_plink <- function(y,
                    file,
                    adjustment = NULL,
                    pop_struct = NULL,
                    verbose = TRUE) {
  ## Made changes to genio::read_bed for our method
  if (missing(file))
    stop('Output file path is required!')
  file <- sub("\\.bed$", "", file)
  require_files_plink(file)
  file_bim = paste0(file, ".bim")
  file_fam = paste0(file, ".fam")
  file_bed = paste0(file, ".bed")

  # Get number of snps/ind
  bim <- genio::read_bim(file_bim, verbose = verbose)
  fam <- genio::read_fam(file_fam, verbose = verbose)
  n_ind <- nrow(fam);
  m_loci <- nrow(bim);

  if (verbose)
    message('Reading: ', file_bed)

  file <- path.expand(file_bed)
  if ( !file.exists( file_bed ) )
    stop( 'File does not exist: ', file_bed )


  if (!(is.matrix(y))) {
    stop("y input must be matrix.")
  }
  if (!(nrow(y) == n_ind)) {
    stop("Observations are rows for y.")
  }
  if (anyNA(y)) {
    stop("Current implementation requires no NAs in trait matrix. Need to adjust inputs accordingly ... feature will be added soon.")
  }

  if (is.null(pop_struct)) {
    h = matrix(1, nrow = nrow(y), ncol = 1)
  } else {
    if (!(is.matrix(pop_struct) & (nrow(y) == nrow(pop_struct)))){
      stop("pop_struct must be matrix and rows are observations.")
    }
    h = cbind(1, pop_struct)
  }

  if (is.null(adjustment)) {
     y <- .quick_lm_cpp(Xs = h, Ys = y)
    } else {
    if (!(is.matrix(adjustment) & (nrow(y) == nrow(adjustment)))){
      stop("pop_struct must be matrix and rows are observations.")
    }
    y <- .quick_lm_cpp(Xs = cbind(h, adjustment), Ys = y)
  }
  lit_out <- .lit_bed_cpp(file = file_bed,
                          m_loci = m_loci,
                          n_ind = n_ind,
                          sY = y,
                          sH = h,
                          verbose = verbose)
  P <- lit_out$P
  colnames(P) <- c("wlit", "ulit", "alit")
  P <- as.data.frame(P);
  bim <- bim[,-3] # drop posg
  bim$maf <- lit_out$MAF
  id <- bim$maf > 0.5
  bim$maf[id] <- 1 - bim$maf[id]
  out <- cbind(bim, P)
  return(out)
}

#' @import stats
#' @import genio
#' @import CompQuadForm
#' @importFrom Rcpp sourceCpp
#' @useDynLib lit
NULL
