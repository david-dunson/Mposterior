// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define SMALL 1e-16

using namespace std;
using namespace Rcpp;
using namespace arma;

/* 
   Calculate log RBF kernel matrix for input matrices x and y. 

   Inputs:
   =======
       x     matrix of size n_x by p_x
       y     matrix of size n_y by p_y
       sigma bandwidth parameter for calculating kernel matrix
   Output:
   =======
       kernmat matrix of dimension n_x by n_y

   Notes:
   ======
   n_x and n_y can be different but p_x must equal p_y.
 */
mat logRbfKernelMatrix (double sigma, mat x, mat y) {
  unsigned int xsize = x.n_rows;
  unsigned int ysize = y.n_rows;

  mat xprod(xsize, xsize), yprod(ysize, ysize), xyprod(xsize, ysize), kernmat(xsize, ysize); 
  xprod.fill(0.0); yprod.fill(0.0); xyprod.fill(0.0); kernmat.fill(0.0);
    
  vec vxprod = diagvec(x * trans(x));
  mat(vxprod).set_size(xsize, 1);
  xprod      = repmat(vxprod, 1, ysize);

  vec vyprod = diagvec(y * trans(y));  
  mat(vyprod).set_size(ysize, 1);
  yprod      = repmat(vyprod, 1, xsize);

  xyprod     = x * trans(y);

  kernmat    = -sigma * xprod - sigma * trans(yprod) + 2 * sigma * xyprod;

  return kernmat;
}

/* 
   Calculate the squared-distance of subset posteriors from M-posterior. 

   Inputs:
   =======
       subsetAtomsList list of length n_s and has subset posterior atoms
       medianAtoms     matrix of size n_m by p_m and has M-posterior atoms 
       subsetProbsList list of length n_s and has probabilites of subset posterior atoms
       medianProbs     matrix of size n_m by 1 and has probabilities of M-posterior atoms 
       sigma           bandwidth parameter for calculating kernel matrix
   Output:
   =======
       norms vector of length n_s and has squared-distance of subset posteriors from M-posterior

   Notes:
   ======
   Each subset posterior is a matrix of size n_a by p_a for a = 1 ... n_s, with rows representing 
   atoms and columns representing dimensions.  All subset posteriors and M-posterior have the same
   dimensions p_1 = ... = p_{n_s} = p_m. The number of atoms n_a can differ in different subset posteriors.
 */
colvec updateWeiszfeldWts (double sigma, field<mat> subsetAtomsList, mat medianAtoms, field<mat> subsetProbsList, mat medianProbs) {
  unsigned int ii;
  unsigned int nsubset = subsetAtomsList.n_elem; // # subset posteriors
  unsigned int natom = medianAtoms.n_rows;       // # M-posterior atoms

  colvec norms(nsubset); norms.fill(0.0); // sq-distance between M-posterior and subset posteriors

  double omax;  

  // matrices for storing 1) distances between M-posterior and subset posterior atoms and 2) probabilities of M-posterior and subset posterior atoms
  mat kernmat1, kernmat2(natom, natom), kernmat12, wtmat1, wtmat2(natom, natom), wtmat12, subWt1, medNorm2, subMed12, subsetAtoms, subsetProbs, subsetProbsMat; 

  kernmat2 = logRbfKernelMatrix(sigma, medianAtoms, medianAtoms); 
  wtmat2   = log(medianProbs * trans(medianProbs));
  medNorm2 = kernmat2 + wtmat2;

  // one round of Weiszfeld updates
  for (ii = 0; ii < nsubset; ++ii) {    
    subsetAtoms = subsetAtomsList(ii);
    subsetProbs = subsetProbsList(ii);

    kernmat1  = logRbfKernelMatrix(sigma, subsetAtoms, subsetAtoms);
    wtmat1    = log(subsetProbs * trans(subsetProbs));

    kernmat12 = logRbfKernelMatrix(sigma, subsetAtoms, medianAtoms);
    wtmat12   =  log(subsetProbs * trans(medianProbs));

    subWt1 = kernmat1 + wtmat1;
    subMed12 = kernmat12 + wtmat12;
    
    omax = MAX(MAX(subWt1.max(), medNorm2.max()), subMed12.max());
    norms(ii) = exp(omax) * (accu(exp(subWt1 - omax)) + accu(exp(medNorm2 - omax)) - 2 * accu(exp(subMed12 - omax)));
    if (norms(ii) < SMALL) {
      Rcpp::Rcout << "norm " << (ii + 1) << " is small; truncating to 1e-16." << std::endl;      
      norms(ii) = SMALL;
    }

    subsetAtoms.reset(); subsetProbs.reset(); subsetProbsMat.reset();
    kernmat1.reset(); wtmat1.reset();
    kernmat12.reset(); wtmat12.reset();
  }

  return norms; 
}

/* 
   Find M-posterior given the subset posteriors using Weiszfeld algorithm. 

   Inputs:
   =======
       subsetPosteriorSamplesList list of length n_s and has samples from subset posteriors
       sigma                      bandwidth parameter for calculating kernel matrix
       maxit                      maximum iterations of Weiszfeld algorithm. Default is 50.
   Output:
   =======
       result list of length 4 that records Weiszfeld weights of subset posteriors, 
              their atoms, and the final estimate of Weiszfeld weights.

   Notes:
   ======
   Each subset posterior is a matrix of size n_a by p_a for a = 1 ... n_s, with rows representing 
   atoms and columns representing dimensions.  All subset posteriors have the same dimensions
   p_1 = ... = p_{n_s}. The number of atoms n_a can differ in different subset posteriors.
 */
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List findWeiszfeldMedian (Rcpp::List subsetAtomsList, double sigma, unsigned int maxit, double tol) {
  unsigned int nsubset = subsetAtomsList.size(); // # subset posteriors
  unsigned int ii, jj, kk, ndim, idx;

  field<mat> subsetPosteriorSamplesList(nsubset);
  mat subsetAtoms;
  SEXP subAtoms;
  for (ii = 0; ii < nsubset; ++ii) {
    subAtoms = subsetAtomsList[ii];
    subsetPosteriorSamplesList(ii) = Rcpp::as<arma::mat>(subAtoms);
  }

  ndim = subsetPosteriorSamplesList(0).n_cols; // # dimensions

  colvec snorms(nsubset), norms(nsubset), normsOld(nsubset); // vectors for storing sq-dist between M-posterior and subset posteriors
  colvec natoms(nsubset), weiszfeldWts(nsubset);  // vectors for storing 1) # atoms in and 2) Weiszfeld wts of subset posteriors
  snorms.fill(0.0); norms.fill(0.0); normsOld.fill(0.0); natoms.fill(0.0); weiszfeldWts.fill(0.0);
    
  field<mat> empiricalMeasureProbList(nsubset); // store the prob of atoms in each subset posterior 

  // initialize empirical subset posterior measure 
  for (ii = 0; ii < nsubset; ++ii) {
    empiricalMeasureProbList(ii) = ones(subsetPosteriorSamplesList(ii).n_rows, 1) / subsetPosteriorSamplesList(ii).n_rows;
    natoms(ii) = subsetPosteriorSamplesList(ii).n_rows;
  }

  mat histWts(nsubset, maxit); histWts.fill(0.0);

  // initialize M-posterior atoms and their measures
  mat medianEmpiricalMeasureAtoms = zeros(sum(natoms), ndim);
  mat medianEmpiricalMeasureProbs = zeros(sum(natoms), 1);
  idx = 0;
  for (ii = 0; ii < nsubset; ++ii) {
    for (jj = 0; jj < natoms(ii); ++jj) {
      medianEmpiricalMeasureProbs(idx, 0) = empiricalMeasureProbList(ii)(jj, 0);
      for (kk = 0; kk < ndim; ++kk) {
        medianEmpiricalMeasureAtoms(idx, kk) = subsetPosteriorSamplesList(ii)(jj, kk);
      }
      idx = idx + 1;
    }
  }
  medianEmpiricalMeasureProbs = medianEmpiricalMeasureProbs / accu(medianEmpiricalMeasureProbs);

  for (jj = 0; jj < maxit; ++jj) {
    if ((jj + 1) % 10 == 0) Rcpp::Rcout << "Weiszfeld iteration: " << jj + 1 << std::endl;

    norms = updateWeiszfeldWts(sigma, subsetPosteriorSamplesList, medianEmpiricalMeasureAtoms, empiricalMeasureProbList, medianEmpiricalMeasureProbs);
    snorms = sqrt(norms);
    weiszfeldWts = (ones(nsubset) / snorms) / sum(ones(nsubset) / snorms);

    idx = 0;
    for (ii = 0; ii < nsubset; ++ii) {
      histWts(ii, jj) =  weiszfeldWts(ii);
      for (kk = 0; kk < natoms(ii); ++kk) {
        medianEmpiricalMeasureProbs(idx, 0) = weiszfeldWts(ii) / natoms(ii);
        idx = idx + 1;
      }
    }
    
    if ((sum(abs(norms - normsOld)) / nsubset) < tol & jj > 10) {
      Rcpp::Rcout << "Weiszfeld algorithm coverged at iteration: " << jj + 1 << std::endl;
      break;
    }
    normsOld = norms;
  }  
  
  if (jj < maxit) histWts = histWts.cols(0, jj);

  Rcpp::List result = List::create(Named("natoms", natoms),
                                   Named("weiszfeldWts", weiszfeldWts),
                                   Named("historyWeiszfeldWts", histWts),
                                   Named("medianAtoms", medianEmpiricalMeasureAtoms)
                                   );  
    
  return result; 
}




