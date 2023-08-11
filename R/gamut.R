#' GAMuT
#'
#' @description
#' The \code{GAMuT} function is a kernel-based multivariate association test.
#' Note that our code to process plink files builds from the
#' \code{\link{genio}} R package.
#'
#' @param y matrix of traits (n observations by k traits)
#' @param file path to plink files
#' @param adjustment matrix of covariates to adjust traits
#' @param pop_struct matrix of PCs that captures population structure
#' @param verbose If TRUE (default) print progress.
#'
#' @return
#' A data frame of p-values where the columns are the cross products/squared residuals
#' and the rows are SNPs.
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
#' out <- gamut_plink(Y, file = file)
#'
#' @seealso \code{\link{lit_plink}}, \code{\link{marginal_plink}}
#' @export
gamut_plink <- function(y,
                      file,
                      adjustment = NULL,
                      pop_struct = NULL,
                      verbose = TRUE) {
  ## Made changes to genio::read_bed for our method
  # die if things are missing
  if (missing(file))
    stop('Output file path is required!')
  file <- sub("\\.bed$", "", file)
  require_files_plink(file)
  file_bim = paste0(file, ".bim")
  file_fam = paste0(file, ".fam")
  file_bed = paste0(file, ".bed")

  ## Get number of snps/ind
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

  gamut_out <- .gamut_bed_cpp(file = file_bed, m_loci = m_loci, n_ind = n_ind, sY = y, sH = h, verbose = verbose)
  P <- gamut_out$P
  colnames(P) <- c("EV", "EQ_EV", "cauchy")
  P <- as.data.frame(P);
  bim <- bim[,-3] # drop posg
  bim$maf <- gamut_out$MAF
  id <- bim$maf > 0.5
  bim$maf[id] <- 1 - bim$maf[id]
  out <- cbind(bim, P)
  return(out)
}
