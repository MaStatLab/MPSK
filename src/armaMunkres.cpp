#include "armaMunkres.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// The code for Munkres minimal assignment algorithm is taken from:
// http://gallery.rcpp.org/articles/minimal-assignment/

void step_one(unsigned int &step, arma::mat &cost, 
              const unsigned int &N) 
{    
  for (unsigned int r = 0; r < N; ++r) {
    cost.row(r) -= arma::min(cost.row(r));
  }
  step = 2;
}

void step_two (unsigned int &step, const arma::mat &cost,
               arma::umat &indM, arma::ivec &rcov, 
               arma::ivec &ccov, const unsigned int &N)
{                      
  for (unsigned int r = 0; r < N; ++r) {
    for (unsigned int c = 0; c < N; ++c) {
      if (cost.at(r, c) == 0.0 && rcov.at(r) == 0 && ccov.at(c) == 0) {
        indM.at(r, c)  = 1;
        rcov.at(r)     = 1;
        ccov.at(c) = 1;
        break;                                              
        // Only take the first
        // zero in a row and column
      }
    }
  }
  /* for later reuse */
  rcov.fill(0);
  ccov.fill(0);
  step = 3;
}

void step_three(unsigned int &step, const arma::umat &indM,
                arma::ivec &ccov, const unsigned int &N) 
{
  unsigned int colcount = 0;
  for (unsigned int r = 0; r < N; ++r) {
    for (unsigned int c = 0; c < N; ++c) {
      if (indM.at(r, c) == 1) {
        ccov.at(c) = 1;
      }
    }
  }
  for (unsigned int c = 0; c < N; ++c) {
    if (ccov.at(c) == 1) {
      ++colcount;
    }
  }
  if (colcount == N) {
    step = 7;
  } else {
    step = 4;
  }
}

void find_noncovered_zero(int &row, int &col,
                          const arma::mat &cost, const arma::ivec &rcov, 
                          const arma::ivec &ccov, const unsigned int &N)
{
  unsigned int r = 0;
  unsigned int c;
  bool done = false;
  row = -1;
  col = -1;
  while (!done) {
    c = 0;
    while (true) {
      if (cost.at(r, c) == 0.0 && rcov.at(r) == 0 && ccov.at(c) == 0) {
        row = r;
        col = c;
        done = true;
      }
      ++c;
      if (c == N || done) {
        break;
      }
    }
    ++r;
    if (r == N) {
      done = true;
    }
  }
}

bool star_in_row(int &row, const arma::umat &indM,
                 const unsigned int &N) 
{
  bool tmp = false;
  for (unsigned int c = 0; c < N; ++c) {
    if (indM.at(row, c) == 1) {
      tmp = true;
      break;
    }
  }
  return tmp;
}

void find_star_in_row (const int &row, int &col, 
                       const arma::umat &indM, const unsigned int &N) 
{
  col = -1;
  for (unsigned int c = 0; c < N; ++c) {
    if (indM.at(row, c) == 1) {
      col = c;
    }
  }
}

void step_four (unsigned int &step, const arma::mat &cost,
                arma::umat &indM, arma::ivec &rcov, arma::ivec &ccov,
                int &rpath_0, int &cpath_0, const unsigned int &N) 
{
  int row = -1;
  int col = -1;
  bool done = false;
  while(!done) {
    find_noncovered_zero(row, col, cost, rcov,
                         ccov, N);                                 
    
    if (row == -1) {                                
      done = true;
      step = 6;
    } else {
      /* uncovered zero */
      indM(row, col) = 2;                         
      if (star_in_row(row, indM, N)) {                            
        find_star_in_row(row, col, indM, N);    
        /* Cover the row with the starred zero
         * and uncover the column with the starred
         * zero. 
         */
        rcov.at(row) = 1;                         
        ccov.at(col) = 0;                          
      } else {
        /* No starred zero in row with 
         * uncovered zero 
         */
        done = true;
        step = 5;
        rpath_0 = row;
        cpath_0 = col;
      }            
    }
  }
}

void find_star_in_col (const int &col, int &row,
                       const arma::umat &indM, const unsigned int &N)
{
  row = -1;
  for (unsigned int r = 0; r < N; ++r) {
    if (indM.at(r, col) == 1) {
      row = r;
    }
  }
}

void find_prime_in_row (const int &row, int &col,
                        const arma::umat &indM, const unsigned int &N)
{
  for (unsigned int c = 0; c < N; ++c) {
    if (indM.at(row, c) == 2) {
      col = c;
    }
  }
}

void augment_path (const int &path_count, arma::umat &indM,
                   const arma::imat &path)
{
  for (unsigned int p = 0; p < path_count; ++p) {
    if (indM.at(path(p, 0), path(p, 1)) == 1) {
      indM.at(path(p, 0), path(p, 1)) = 0;
    } else {
      indM.at(path(p, 0), path(p, 1)) = 1;
    }
  }
}

void clear_covers (arma::ivec &rcov, arma::ivec &ccov)
{
  rcov.fill(0);
  ccov.fill(0);
}

void erase_primes(arma::umat &indM, const unsigned int &N)
{
  for (unsigned int r = 0; r < N; ++r) {
    for (unsigned int c = 0; c < N; ++c) {
      if (indM.at(r, c) == 2) {
        indM.at(r, c) = 0;
      }
    }
  }
}

void step_five (unsigned int &step,
                arma::umat &indM, arma::ivec &rcov, 
                arma::ivec &ccov, arma::imat &path, 
                int &rpath_0, int &cpath_0, 
                const unsigned int &N)
{
  bool done = false;
  int row = -1;
  int col = -1;
  unsigned int path_count = 1;
  path.at(path_count - 1, 0) = rpath_0;
  path.at(path_count - 1, 1) = cpath_0;
  while (!done) {
    find_star_in_col(path.at(path_count - 1, 1), row, 
                     indM, N);
    if (row > -1) {                                
      /* Starred zero in row 'row' */
      ++path_count;
      path.at(path_count - 1, 0) = row;
      path.at(path_count - 1, 1) = path.at(path_count - 2, 1);
    } else {
      done = true;
    }
    if (!done) {
      /* If there is a starred zero find a primed 
       * zero in this row; write index to 'col' */
      find_prime_in_row(path.at(path_count - 1, 0), col, 
                        indM, N);  
      ++path_count;
      path.at(path_count - 1, 0) = path.at(path_count - 2, 0);
      path.at(path_count - 1, 1) = col;
    }
  }
  augment_path(path_count, indM, path);
  clear_covers(rcov, ccov);
  erase_primes(indM, N);
  step = 3;
}

void find_smallest (double &minval, const arma::mat &cost, 
                    const arma::ivec &rcov, const arma::ivec &ccov, 
                    const unsigned int &N)
{
  for (unsigned int r = 0; r < N; ++r) {
    for (unsigned int c = 0; c < N; ++c) {
      if (rcov.at(r) == 0 && ccov.at(c) == 0) {                                                                    
        if (minval > cost.at(r, c)) {
          minval = cost.at(r, c);
        }
      }
    }
  }
}

void step_six (unsigned int &step, arma::mat &cost,
               const arma::ivec &rcov, const arma::ivec &ccov, 
               const unsigned int &N) 
{
  // double minval = DBL_MAX;
  double minval = DOUBLE_XMAX;
  find_smallest(minval, cost, rcov, ccov, N);
  for (unsigned int r = 0; r < N; ++r) {
    for (unsigned int c = 0; c < N; ++c) {
      if (rcov.at(r) == 1) {
        cost.at(r, c) += minval;
      }
      if (ccov.at(c) == 0) {
        cost.at(r, c) -= minval;
      }
    }
  }
  step = 4;
}

arma::umat hungarian(const arma::mat &input_cost)
{
  const unsigned int N = input_cost.n_rows;
  unsigned int step = 1;
  int cpath_0 = 0;
  int rpath_0 = 0;
  arma::mat cost(input_cost);
  arma::umat indM(N, N);
  arma::ivec rcov(N);
  arma::ivec ccov(N);
  arma::imat path(2 * N, 2);
  
  indM = arma::zeros<arma::umat>(N, N);
  bool done = false;
  while (!done) {
    switch (step) {
    case 1:
      step_one(step, cost, N);
      break;
    case 2:
      step_two(step, cost, indM, rcov, ccov, N);
      break;
    case 3:
      step_three(step, indM, ccov, N);
      break;
    case 4:
      step_four(step, cost, indM, rcov, ccov,
                rpath_0, cpath_0, N);
      break;
    case 5: 
      step_five(step, indM, rcov, ccov,
                path, rpath_0, cpath_0, N);
      break;
    case 6:
      step_six(step, cost, rcov, ccov, N);
      break;
    case 7:            
      done = true;
      break;
    }
  }
  return indM;
}

// [[Rcpp::export]]
arma::umat hungarian_cc(arma::mat cost)
{
  arma::umat indM = hungarian(cost);
  return(indM);
}

// [[Rcpp::export]]
Rcpp::List relabel(const Rcpp::List res)
{
  Rcpp::List chain = clone(Rcpp::as<Rcpp::List>(res["chain"]));
  Rcpp::List prior = Rcpp::as<Rcpp::List>(res["prior"]);
  int K = Rcpp::as<int>(prior["K"]);
  cout << "K=" << K << endl;
  umat t = Rcpp::as<umat>(chain["t"]);
  int T = t.n_rows;
  cout << "T=" << T << endl;
  int N = t.n_cols;
  cout << "N=" << N << endl;
  urowvec refZ = t.row(T-1);

  umat t_relabel = Rcpp::as<umat>(chain["t"]);
  // cout << unique(t_relabel) << endl;
  // mat z_relabel = Rcpp::as<mat>(chain["z"]);
  cube W_relabel = Rcpp::as<cube>(chain["W"]);
  int J = W_relabel.n_rows;
  cout << "J=" << J << endl;
  cube xi_relabel = Rcpp::as<cube>(chain["xi"]);
  int p = xi_relabel.n_cols / K;
  cout << "p=" << p << endl;
  cube xi0_relabel = Rcpp::as<cube>(chain["xi0"]);
  cube psi_relabel = Rcpp::as<cube>(chain["psi"]);
  cube G_relabel = Rcpp::as<cube>(chain["G"]);
  cube E_relabel = Rcpp::as<cube>(chain["E"]);
  umat S_relabel = Rcpp::as<umat>(chain["S"]);
  mat log_py_relabel = Rcpp::as<mat>(chain["log_py"]);
  mat perplexity_relabel = Rcpp::as<mat>(chain["perplexity"]);
  mat nResampled_relabel = Rcpp::as<mat>(chain["nResampled"]);
  cube Omega_relabel = Rcpp::as<cube>(chain["Omega"]);
  cube alpha_relabel = Rcpp::as<cube>(chain["alpha"]);

  uvec refZobs(K);
  refZobs.fill(0);
  for (int i=0; i<N; i++) {
    refZobs(refZ(i))++;
  }

  mat cost(K,K);
  umat permut;
  
  for (int s=0; s<(T-1); s++) {
  // for (int s=0; s<1; s++) {
    urowvec currZ = t.row(s);
    int refClass, curClass;
    // Find relabeling cost matrix:
    for (int h=0; h<K; h++) {
      curClass = refZobs(h);
      for (int j=0; j<K; j++) {
        cost(h,j) = curClass;
      }
    }
  
    for (int i=0; i<N; i++) {
      refClass = refZ(i);
      curClass = currZ(i);
      cost(refClass,curClass) -= 1.0;
    }
  
    permut = hungarian_cc(cost);
    // cout << permut << endl;

    urowvec relabeling(K);
    // urowvec back_relabeling(K);
    for (int row=0; row<K; row++) {
      for (int col=0; col<K; col++) {
        if (permut(row,col)==1) {
          relabeling(row) = col;
          // back_relabeling(col) = row;
        }
      }
    }
    
    // cout << relabeling << endl;
    // cout << relabeling(fined(relabeling==t_copy(i))) << endl;
    
    
    urowvec t_copy = t_relabel.row(s);
    for (int i=0; i<N; i++) {
      // t_relabel(s, i) = back_relabeling(t_copy(i));
      uvec idx = find(relabeling==t_copy(i));
      t_relabel(s, i) = relabeling(idx(0));
    }
    

    mat W_copy = W_relabel.slice(s);
    mat xi_copy = xi_relabel.slice(s);
    mat xi0_copy = xi0_relabel.slice(s);
    mat psi_copy = psi_relabel.slice(s);
    mat G_copy = G_relabel.slice(s);
    mat E_copy = E_relabel.slice(s);
    urowvec S_copy = S_relabel.row(s);
    rowvec log_py_copy = log_py_relabel.row(s);
    rowvec perplexity_copy = perplexity_relabel.row(s);
    rowvec nResampled_copy = nResampled_relabel.row(s);
    mat Omega_copy = Omega_relabel.slice(s);
    mat alpha_copy = alpha_relabel.slice(s);
    for (int k=0; k<K; k++) {
      uvec idx = find(relabeling==k);
      int kk = relabeling(idx(0));
      W_relabel.slice(s).col(k) = W_copy.col(kk);
      xi_relabel.subcube(0,k*p,s,J-1,k*p+p-1,s) = 
        xi_copy.submat(0,kk*p,J-1,kk*p+p-1);
      xi0_relabel.slice(s).col(k) = xi0_copy.col(kk);
      psi_relabel.slice(s).col(k) = psi_copy.col(kk);
      G_relabel.subcube(0,k*p,s,p-1,k*p+p-1,s) = 
        G_copy.submat(0,kk*p,p-1,kk*p+p-1);
      E_relabel.subcube(0,k*p,s,p-1,k*p+p-1,s) = 
        E_copy.submat(0,kk*p,p-1,kk*p+p-1);
      // S_relabel.slice(s).col(k) = S_copy.col(kk);
      S_relabel(s, k) = S_copy(kk);
      log_py_relabel(s, k) = log_py_copy(kk);
      perplexity_relabel(s, k) = perplexity_copy(kk);
      nResampled_relabel(s, k) = nResampled_copy(kk);
      Omega_relabel.subcube(0,k*p,s,p-1,k*p+p-1,s) = 
        Omega_copy.submat(0,kk*p,p-1,kk*p+p-1);
      alpha_relabel.slice(s).col(k) = alpha_copy.col(kk);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named( "t" ) = t_relabel,
    Rcpp::Named( "W" ) = W_relabel,
    Rcpp::Named( "xi" ) = xi_relabel,
    Rcpp::Named( "xi0" ) = xi0_relabel,
    Rcpp::Named( "psi" ) = psi_relabel,
    Rcpp::Named( "G" ) = G_relabel,
    Rcpp::Named( "E" ) = E_relabel,
    Rcpp::Named( "S" ) = S_relabel,
    Rcpp::Named( "varphi" ) = Rcpp::as<vec>(chain["varphi"]),
    Rcpp::Named( "a0" ) = Rcpp::as<vec>(chain["a0"]),
    Rcpp::Named( "log_py" ) = log_py_relabel,
    Rcpp::Named( "perplexity" ) = perplexity_relabel,
    Rcpp::Named( "nResampled" ) = nResampled_relabel,
    Rcpp::Named( "Omega" ) = Omega_relabel,
    Rcpp::Named( "alpha" ) = alpha_relabel
  );
}
