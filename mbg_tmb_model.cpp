// ///////////////////////////////////////////////////
// NOTES:
// 1. Type `Type` is a special TMB type that should be used for all variables, except `int` can be used as well
// 2. Anything in the density namespace (ie SEPARABLE) returns the negative log likelihood and is thus added to accumulator
//    also, other density function such as dnorm and dbinom return positive log likelihood and are thus subtracted away.
// 3. ref https://github.com/nmmarquez/re_simulations/blob/master/inla/sta.cpp
//        https://github.com/nmmarquez/re_simulations/blob/master/inla/SPDEAR1AR1.R
// ///////////////////////////////////////////////////

// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
    SparseMatrix<Type> Q;
    Type kappa2 = exp(2. * logkappa);
    Type kappa4 = kappa2*kappa2;
    Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
    return Q;
}

// Robust Inverse Logit that sets min and max values to avoid numerical instability
template<class Type>
Type invlogit_robust(Type x){
  if (x < -20.723){
    x = -20.723; // corresponds to p=1e-9
  } else if ( x > 20.723 ){
    x = 20.723;  // cooresponds to p=1-1e-9
  }
  return 1 / (1 + exp( -1.0 * x ));
}


// Corresponding list object on the C++ side
template<class Type>
struct option_list {
  int use_priors;
  int adreport_off;
  int use_poisson;
  // Way easier to read these in as vectors of integers and then index the first
  option_list(SEXP x){
    use_priors = asVector<int>(getListElement(x,"use_priors"))[0];
    adreport_off = asVector<int>(getListElement(x,"adreport_off"))[0];
    use_poisson = asVector<int>(getListElement(x,"use_poisson"))[0];
  }
};


// Constrain alpha values
template<class Type>
vector<Type> constrain_pars(vector<Type> alpha, vector<int> constraints){
  int K = alpha.size();
  vector<Type> alpha_c(K);
  
  for(int k = 0; k < K; k++){
    if(constraints[k] == 1){
      alpha_c[k] = exp(alpha[k]);
    }
    if(constraints[k] == -1){
      alpha_c[k] = -1. * exp(alpha[k]);
    }
    if(constraints[k] == 0){
      alpha_c[k] = alpha[k];
    }
  }
  
  return alpha_c;
}


// objective function (ie the likelihood function for the model), returns the evaluated negative log likelihood
template<class Type>
Type objective_function<Type>::operator() ()
{

  // ////////////////////////////////////////////////////////////////////////////
  // INPUTS
  // ////////////////////////////////////////////////////////////////////////////
  DATA_INTEGER(flag); // flag=0 => only prior

  // Indices
  DATA_INTEGER(num_i);       // number of datapts
  DATA_INTEGER(num_s);       // number of mesh pts in space mesh

  // Data (each, excpect for X_ij is a vector of length num_i)
  DATA_VECTOR(y_i);          // obs successes per binomial experiment at point i (aka cluster)
  DATA_VECTOR(n_i);          // trials per cluster
  DATA_IVECTOR(w_i);         // weights for observations
  DATA_MATRIX(X_ij);         // covariate design matrix (num_i by number of fixed effects matrix)
  DATA_IVECTOR(fconstraints); // constraints of fixed effects

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M1);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M2);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(Aproj); // used to project from spatial mesh to data locations

  // Options
  DATA_STRUCT(options, option_list);  // boolean vector of options to be used to select different models/modelling options: see above
  
  // Parameters
  PARAMETER_VECTOR(alpha_j);   // fixed effect coefs, including intercept as first index
  PARAMETER(logtau);           // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(logkappa);         // log of INLA kappa - related to spatial correlation and range

  // Random effects
  PARAMETER_VECTOR(epsilon_s); // Random effects for each space mesh location.

  printf("Epsilon_s size: %ld \n", epsilon_s.size());

  // ////////////////////////////////////////////////////////////////////////////
  // LIKELIHOOD
  // ////////////////////////////////////////////////////////////////////////////

  // Define the joint-negative log-likelihood as a parallel_accumulator
  // this allows us to add or subtract numbers to the object in parallel
  // parallel_accumulator<Type> jnll(this);
  Type jnll = 0;

  // print parallel info
  max_parallel_regions = omp_get_max_threads();

  // Make spatial, t and z precision matrices
  SparseMatrix<Type> Q_ss   = spde_Q(logkappa, logtau, M0, M1, M2);
  printf("Q_ss size: %ld \n", Q_ss.size());

  // Make transformations of some of our parameters
  Type range         = sqrt(8.0) / exp(logkappa);
  Type sigma         = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) * exp(2.0 * logkappa));

  // Define objects for derived values
  vector<Type> fe_i(num_i);                         // main effect X_ij %*% t(alpha_j)
  vector<Type> projepsilon_i(num_i);                // value of gmrf at data points
  vector<Type> prob_i(num_i);                       // Logit estimated prob for each point i

  // Prior contribution to likelihood. Values are defaulted.
  // Only run if options.use_priors==1
  if(options.use_priors == 1) {
   PARALLEL_REGION jnll -= dnorm(logtau,    Type(0.0), Type(1.0), true);  // N(0,1) prior for logtau
   PARALLEL_REGION jnll -= dnorm(logkappa,  Type(0.0), Type(1.0), true);  // N(0,1) prior for logkappa
   // Additional constraint on sigma
   PARALLEL_REGION jnll -= dnorm(sigma,     Type(0.0), Type(0.25), true);
   for( int j = 0; j < alpha_j.size(); j++){
     PARALLEL_REGION jnll -= dnorm(alpha_j(j), Type(0.0), Type(3), true); // N(0, sqrt(1/.001)) prior for fixed effects.
   }
  }

  // Latent field/Random effect contribution to likelihood.
  printf("GP FOR SPACE  ONLY \n");
  PARALLEL_REGION jnll += GMRF(Q_ss,false)(epsilon_s);

  // Project from mesh points to data points in order to eval likelihood at each data point
  printf("Project Epsilon \n");
  projepsilon_i = Aproj * epsilon_s.matrix();

  // evaluate fixed effects for alpha_j values
  vector<Type> calpha_j = constrain_pars(alpha_j, fconstraints);
  fe_i = X_ij * calpha_j.matrix();

  // Return un-normalized density on request
  if (flag == 0) return jnll;

  // Likelihood contribution from each datapoint i
  printf("Data likelihood \n");
  for (int i = 0; i < num_i; i++){
    // mean model
    prob_i(i) = fe_i(i) + projepsilon_i(i);
    // Optionally switch between a Poisson and a Binomial likelihood
    if( options.use_poisson == 1 ) {
      // Poisson likelihood with a log link
      PARALLEL_REGION jnll -= dpois( y_i(i), n_i(i) * exp(prob_i(i)), true ) * w_i(i);
    } else {
      // Binomia likelihood with a logit link
      PARALLEL_REGION jnll -= dbinom( y_i(i), n_i(i), invlogit_robust(prob_i(i)), true ) * w_i(i);      
    }
  }

  // Report estimates
  if(options.adreport_off == 0){
    ADREPORT(alpha_j);
  }

  return jnll;
}
