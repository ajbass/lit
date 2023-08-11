// The code to input SNPs from a bed file is from Alex Ochoa's genio package,
// where changes were made to integrate the LIT method.
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cerrno>
#include "lit_cpp.h"
#include "marginal.h"
#include "linreg.h"
#include "svd.h"
#include "hsic_cpp.h"
#include "gamut.h"

using namespace Rcpp;

// expected header (magic numbers)
// assume standard locus-major order and latest format
const unsigned char plink_bed_byte_header[3] = {0x6c, 0x1b, 1};

void flush_console() {
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
  R_FlushConsole();
#endif
}

/**
 * Latent Interaction Testing
 *
 * @param file location of bed file
 * @param m_loci total number of loci
 * @param n_ind sample size
 * @param sY traits
 * @param sH population structure
 * @param tot total number of SQ/CP terms
 * @param verbose print progress
 *
 * @return p-values
 */
// [[Rcpp::export(.lit_bed_cpp)]]
List lit_bed_cpp(const char* file, int m_loci, int n_ind, arma::mat sY,
                 arma::mat sH, bool verbose) {
  // - file assumed to be full path (no missing extensions)
  // unfortunately BED format requires dimensions to be known
  // (so outside this function, the BIM and FAM files must be parsed first)
  // Obtain environment containing function
  // open input file in "binary" mode
  FILE *file_stream = fopen( file, "rb" );
  // die right away if needed, before initializing buffers etc
  if ( file_stream == NULL ) {
    // send error message to R
    //    char msg[100];
    Rprintf("Truncated file"); // convert to 1-based coordinates
    //      stop(msg);
  }
  /////////////////////////
  // check magic numbers //
  /////////////////////////
  double indexTmp =0;
 // std::cout << "[";
  // for header only
  unsigned char *buffer_header = (unsigned char *) malloc( 3 );

  // for extra sanity checks, keep track of bytes actually read (to recognize truncated files)
  // reuse this one for genotypes below
  size_t n_buf_read;

  // read header bytes (magic numbers)
  n_buf_read = fread( buffer_header, sizeof(unsigned char), 3, file_stream );
  // this might just indicate an empty file
  if ( n_buf_read != 3 ) {
    // wrap up everything properly
    free( buffer_header ); // free buffer memory
    fclose( file_stream ); // close file
    // now send error message to R
    stop("Input BED file did not have a complete header (3-byte magic numbers)!");
  }

  // require that they match our only supported specification of locus-major order and latest format
  // was using strcmp but there are funky issues (wants signed, but we don't really want order anyway, just test for equality)
  // use explicit loop instead
  int pos;
  for (pos = 0; pos < 3; pos++) {
    if ( plink_bed_byte_header[pos] != buffer_header[pos] ) {
      // wrap up everything properly
      free( buffer_header ); // free buffer memory
      fclose( file_stream ); // close file
      // now send error message to R
      stop("Input BED file is not in supported format.  Either magic numbers do not match, or requested sample-major format is not supported.  Only latest locus-major format is supported!");
    }
  }

  // free header buffer, completely done with it
  free( buffer_header );

  ////////////////////
  // read genotypes //
  ////////////////////

  // number of columns (bytes) in input (for buffer), after byte compression
  // size set for full row, but overloaded used first for this header comparison
  // chose size_t to have it match n_buf_read value returned by fread
  size_t n_buf = ( n_ind + 3 ) / 4;
  // initialize row buffer
  unsigned char *buffer = (unsigned char *) malloc( n_buf );
  // navigate data and process
  NumericMatrix P (m_loci, 3);
  NumericVector MAF (m_loci);
  int i;
  size_t k; // to match n_buf type
  unsigned char buf_k; // working of buffer at k'th position
  unsigned char xij; // copy of extracted genotype

  arma::vec X (n_ind, arma::fill::zeros);
  arma::colvec index (n_ind, arma::fill::ones);
  for (i = 0; i < m_loci; i++) {
    // read whole row into buffer
    n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );

    // always check that file was not done too early
    if ( n_buf_read != n_buf ) {
      // wrap up everything properly
      free( buffer ); // free buffer memory
      fclose( file_stream ); // close file
      // now send error message to R
      //    char msg[100];
      Rprintf("Truncated file"); // convert to 1-based coordinates
      //      stop(msg);
    }

    // process buffer now!

    // always reset these at start of row
    int j = 0; // individuals

    // navigate buffer positions k (not individuals j)
    for (k = 0; k < n_buf; k++) {
      // copy down this value, which will be getting edited
      buf_k = buffer[k];

      // navigate the four positions
      // pos is just a dummy counter not really used except to know when to stop
      // update j too, accordingly
      for (pos = 0; pos < 4; pos++, j++) {

        if (j < n_ind) {
          // extract current genotype using this mask
          // (3 == 00000011 in binary)
          xij = buf_k & 3;

          // re-encode into proper values, store in R matrix
          // for maximum speed, test for most common cases first: the homozygotes
          // - this is because of the binomial expansion:
          //   (p-q)^2 = p^2 + q^2 - 2pq > 0,
          //   so `2pq` is always smaller than `p^2 + q^2`.
          //   `2pq` becomes rarer under population structure!
          // next most common ought to be the heterozygote, then NA, but as these are mutually exclusive then it doesn't matter
          if (xij == 0) {
            X(j) = 2; // 0 -> 2
          } else if (xij == 2) {
            X(j) = 1; // 2 -> 1
          } else if (xij == 1) {
            index(j) = 0;
          }
          // there are no other values, so 3 -> 0 must be the case here
          // R's IntegerMatrix are initialized to zeroes, so no edits are necessary!

          // shift packed data, throwing away genotype we just processed
          buf_k = buf_k >> 2;
        }
      }
    }
    arma::uvec index2 = find(index);
    arma::vec X2 = X.rows(index2);
    arma::mat H2 = sH.rows(index2);
    arma::mat Y2 = sY.rows(index2);
    MAF(i) = mean(X2) / 2;
    P( i , _ ) = lit_internal(X2, Y2, H2);
    index.ones();
    if (i > indexTmp && verbose) {
      Rcpp::Rcout << int(indexTmp / m_loci * 100.0) << "%\r";
      flush_console();
      indexTmp = indexTmp + m_loci / 100.0;
    }
    X.zeros();
  }
  // finished matrix!
  // let's check that file was indeed done
  n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );
  // wrap up regardless
  // and more troubleshooting messages (for windows)
  if ( fclose( file_stream ) != 0 )
    stop("Input BED file stream close failed!");
  free( buffer );
  if ( n_buf_read != 0 ) {
    // now send error message to R
    stop("Input BED file continued after all requested rows were read!  Either the specified the number of loci was too low or the input file is corrupt!");
  }
  List ret;
  ret["P"] = P;
  ret["MAF"] = MAF;
  return(ret);
}

/**
 * GAMuT
 *
 * @param file location of bed file
 * @param m_loci total number of loci
 * @param n_ind sample size
 * @param sY traits
 * @param sH population structure
 * @param verbose print progress
 * @return p-values
 */
// [[Rcpp::export(.gamut_bed_cpp)]]
List gamut_bed_cpp(const char* file, int m_loci, int n_ind, arma::mat sY,
                   arma::mat sH, bool verbose) {
  // - file assumed to be full path (no missing extensions)
  // unfortunately BED format requires dimensions to be known
  // (so outside this function, the BIM and FAM files must be parsed first)
  // Obtain environment containing function
  // open input file in "binary" mode
  double indexTmp = 0;
  FILE *file_stream = fopen( file, "rb" );
  // die right away if needed, before initializing buffers etc
  if ( file_stream == NULL ) {
    // send error message to R
    //    char msg[100];
    Rprintf("Truncated file"); // convert to 1-based coordinates
    //      stop(msg);
  }
  /////////////////////////
  // check magic numbers //
  /////////////////////////

  // for header only
  unsigned char *buffer_header = (unsigned char *) malloc( 3 );

  // for extra sanity checks, keep track of bytes actually read (to recognize truncated files)
  // reuse this one for genotypes below
  size_t n_buf_read;

  // read header bytes (magic numbers)
  n_buf_read = fread( buffer_header, sizeof(unsigned char), 3, file_stream );
  // this might just indicate an empty file
  if ( n_buf_read != 3 ) {
    // wrap up everything properly
    free( buffer_header ); // free buffer memory
    fclose( file_stream ); // close file
    // now send error message to R
    stop("Input BED file did not have a complete header (3-byte magic numbers)!");
  }

  // require that they match our only supported specification of locus-major order and latest format
  // was using strcmp but there are funky issues (wants signed, but we don't really want order anyway, just test for equality)
  // use explicit loop instead
  int pos;
  for (pos = 0; pos < 3; pos++) {
    if ( plink_bed_byte_header[pos] != buffer_header[pos] ) {
      // wrap up everything properly
      free( buffer_header ); // free buffer memory
      fclose( file_stream ); // close file
      // now send error message to R
      stop("Input BED file is not in supported format.  Either magic numbers do not match, or requested sample-major format is not supported.  Only latest locus-major format is supported!");
    }
  }

  // free header buffer, completely done with it
  free( buffer_header );


  ////////////////////
  // read genotypes //
  ////////////////////

  // number of columns (bytes) in input (for buffer), after byte compression
  // size set for full row, but overloaded used first for this header comparison
  // chose size_t to have it match n_buf_read value returned by fread
  size_t n_buf = ( n_ind + 3 ) / 4;
  // initialize row buffer
  unsigned char *buffer = (unsigned char *) malloc( n_buf );
  // navigate data and process
  NumericMatrix P (m_loci, 3);
  NumericVector MAF (m_loci);
  int i;
  size_t k; // to match n_buf type
  unsigned char buf_k; // working of buffer at k'th position
  unsigned char xij; // copy of extracted genotype

  arma::vec X (n_ind, arma::fill::zeros);
  arma::colvec index (n_ind, arma::fill::ones);
  for (i = 0; i < m_loci; i++) {
    // read whole row into buffer
    n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );

    // always check that file was not done too early
    if ( n_buf_read != n_buf ) {
      // wrap up everything properly
      free( buffer ); // free buffer memory
      fclose( file_stream ); // close file
      // now send error message to R
  //    char msg[100];
      Rprintf("Truncated file"); // convert to 1-based coordinates
//      stop(msg);
    }

    // process buffer now!

    // always reset these at start of row
    int j = 0; // individuals

    // navigate buffer positions k (not individuals j)
    for (k = 0; k < n_buf; k++) {
      // copy down this value, which will be getting edited
      buf_k = buffer[k];

      // navigate the four positions
      // pos is just a dummy counter not really used except to know when to stop
      // update j too, accordingly
      for (pos = 0; pos < 4; pos++, j++) {

        if (j < n_ind) {
          // extract current genotype using this mask
          // (3 == 00000011 in binary)
          xij = buf_k & 3;

          // re-encode into proper values, store in R matrix
          // for maximum speed, test for most common cases first: the homozygotes
          // - this is because of the binomial expansion:
          //   (p-q)^2 = p^2 + q^2 - 2pq > 0,
          //   so `2pq` is always smaller than `p^2 + q^2`.
          //   `2pq` becomes rarer under population structure!
          // next most common ought to be the heterozygote, then NA, but as these are mutually exclusive then it doesn't matter
          if (xij == 0) {
            X(j) = 2; // 0 -> 2
          } else if (xij == 2) {
            X(j) = 1; // 2 -> 1
          } else if (xij == 1) {
            index(j) = 0;
          }
          // there are no other values, so 3 -> 0 must be the case here
          // R's IntegerMatrix are initialized to zeroes, so no edits are necessary!

          // shift packed data, throwing away genotype we just processed
          buf_k = buf_k >> 2;
        }
      }
    }
    arma::uvec index2 = find(index);
    arma::vec X2 = X.rows(index2);
    arma::mat H2 = sH.rows(index2);
    arma::mat Y2 = sY.rows(index2);
    MAF(i) = mean(X2) / 2;
    P( i , _ ) = gamut_internal(X2, Y2, H2);
    index.ones();
    if (i > indexTmp && verbose) {
      Rcpp::Rcout << int(indexTmp / m_loci * 100.0) << "%\r";
      flush_console();
      indexTmp = indexTmp + m_loci / 100.0;
    }
    X.zeros();
  }
  // finished matrix!
  // let's check that file was indeed done
  n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );
  // wrap up regardless
  // and more troubleshooting messages (for windows)
  if ( fclose( file_stream ) != 0 )
    stop("Input BED file stream close failed!");
  free( buffer );
  if ( n_buf_read != 0 ) {
    // now send error message to R
    stop("Input BED file continued after all requested rows were read!  Either the specified the number of loci was too low or the input file is corrupt!");
  }
  List ret;
  ret["P"] = P;
  ret["MAF"] = MAF;
  return(ret);
}

/**
 * Marginal Testing
 *
 * @param file location of bed file
 * @param m_loci total number of loci
 * @param n_ind sample size
 * @param sY traits
 * @param sH population structure
 * @param verbose print progress
 *
 * @return p-values
 */
// [[Rcpp::export(.marginal_bed_cpp)]]
List marginal_bed_cpp(const char* file, int m_loci, int n_ind, arma::mat sY,
                      arma::mat sH, int tot, bool verbose) {
  // - file assumed to be full path (no missing extensions)
  // unfortunately BED format requires dimensions to be known
  // (so outside this function, the BIM and FAM files must be parsed first)
  // Obtain environment containing function
  // open input file in "binary" mode
  FILE *file_stream = fopen( file, "rb" );
  double indexTmp = 0;
  // die right away if needed, before initializing buffers etc
  if ( file_stream == NULL ) {
    // send error message to R
    //    char msg[100];
    Rprintf("Could not open BED file"); // convert to 1-based coordinates
    //      stop(msg);
  }
  /////////////////////////
  // check magic numbers //
  /////////////////////////

  // for header only
  unsigned char *buffer_header = (unsigned char *) malloc( 3 );

  // for extra sanity checks, keep track of bytes actually read (to recognize truncated files)
  // reuse this one for genotypes below
  size_t n_buf_read;

  // read header bytes (magic numbers)
  n_buf_read = fread( buffer_header, sizeof(unsigned char), 3, file_stream );
  // this might just indicate an empty file
  if ( n_buf_read != 3 ) {
    // wrap up everything properly
    free( buffer_header ); // free buffer memory
    fclose( file_stream ); // close file
    // now send error message to R
    stop("Input BED file did not have a complete header (3-byte magic numbers)!");
  }

  // require that they match our only supported specification of locus-major order and latest format
  // was using strcmp but there are funky issues (wants signed, but we don't really want order anyway, just test for equality)
  // use explicit loop instead
  int pos;
  for (pos = 0; pos < 3; pos++) {
    if ( plink_bed_byte_header[pos] != buffer_header[pos] ) {
      // wrap up everything properly
      free( buffer_header ); // free buffer memory
      fclose( file_stream ); // close file
      // now send error message to R
      stop("Input BED file is not in supported format.  Either magic numbers do not match, or requested sample-major format is not supported.  Only latest locus-major format is supported!");
    }
  }

  // free header buffer, completely done with it
  free( buffer_header );


  ////////////////////
  // read genotypes //
  ////////////////////

  // number of columns (bytes) in input (for buffer), after byte compression
  // size set for full row, but overloaded used first for this header comparison
  // chose size_t to have it match n_buf_read value returned by fread
  size_t n_buf = ( n_ind + 3 ) / 4;
  // initialize row buffer
  unsigned char *buffer = (unsigned char *) malloc( n_buf );
  // navigate data and process
  NumericMatrix P (m_loci, tot);
  NumericVector MAF (m_loci);
  int i;
  size_t k; // to match n_buf type
  unsigned char buf_k; // working of buffer at k'th position
  unsigned char xij; // copy of extracted genotype

  arma::vec X (n_ind, arma::fill::zeros);
  arma::colvec index (n_ind, arma::fill::ones);
  for (i = 0; i < m_loci; i++) {
    // read whole row into buffer
    n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );

    // always check that file was not done too early
    if ( n_buf_read != n_buf ) {
      // wrap up everything properly
      free( buffer ); // free buffer memory
      fclose( file_stream ); // close file
      // now send error message to R
      Rprintf("Truncated file"); // convert to 1-based coordinates
    }

    // process buffer now!

    // always reset these at start of row
    int j = 0; // individuals

    // navigate buffer positions k (not individuals j)
    for (k = 0; k < n_buf; k++) {
      // copy down this value, which will be getting edited
      buf_k = buffer[k];

      // navigate the four positions
      // pos is just a dummy counter not really used except to know when to stop
      // update j too, accordingly
      for (pos = 0; pos < 4; pos++, j++) {

        if (j < n_ind) {
          // extract current genotype using this mask
          // (3 == 00000011 in binary)
          xij = buf_k & 3;

          // re-encode into proper values, store in R matrix
          // for maximum speed, test for most common cases first: the homozygotes
          // - this is because of the binomial expansion:
          //   (p-q)^2 = p^2 + q^2 - 2pq > 0,
          //   so `2pq` is always smaller than `p^2 + q^2`.
          //   `2pq` becomes rarer under population structure!
          // next most common ought to be the heterozygote, then NA, but as these are mutually exclusive then it doesn't matter
          if (xij == 0) {
            X(j) = 2; // 0 -> 2
          } else if (xij == 2) {
            X(j) = 1; // 2 -> 1
          } else if (xij == 1) {
            index(j) = 0;
          }
          // there are no other values, so 3 -> 0 must be the case here
          // R's IntegerMatrix are initialized to zeroes, so no edits are necessary!

          // shift packed data, throwing away genotype we just processed
          buf_k = buf_k >> 2;
        }
      }
    }
    arma::uvec index2 = find(index);
    arma::vec X2 = X.rows(index2);
    arma::mat H2 = sH.rows(index2);
    arma::mat Y2 = sY.rows(index2);
    MAF(i) = mean(X2) / 2;
    P( i , _ ) = marginal_internal(X2, Y2, H2);
    index.ones();
    if (i > indexTmp && verbose) {
      Rcpp::Rcout << int(indexTmp / m_loci * 100.0) << "%\r";
      flush_console();
      indexTmp = indexTmp + m_loci / 100.0;
    }
    X.zeros();
  }
  // finished matrix!
  // let's check that file was indeed done
  n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );
  // wrap up regardless
  // and more troubleshooting messages (for windows)
  if ( fclose( file_stream ) != 0 )
    stop("Input BED file stream close failed!");
  free( buffer );
  if ( n_buf_read != 0 ) {
    // now send error message to R
    stop("Input BED file continued after all requested rows were read!  Either the specified the number of loci was too low or the input file is corrupt!");
  }
  List ret;
  ret["P"] = P;
  ret["MAF"] = MAF;
  return(ret);
}
