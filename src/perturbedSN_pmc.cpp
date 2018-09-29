#include "RcppArmadillo.h"
#include "perturbedSN_pmc.h"
#include "perturbedSN_helpers.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

const double log2pi = std::log(2.0 * M_PI);

PMC::PMC(arma::mat Y,
         arma::uvec C,
         Rcpp::List prior,
         Rcpp::List pmc,
         Rcpp::List state ) : Y(Y), C(C)
{
  p = Y.n_cols;
  n = Y.n_rows;
  J = C.max() + 1;
  K = Rcpp::as<int>(prior["K"]);

  num_particles = Rcpp::as<int>(pmc["npart"]);
  num_iter = Rcpp::as<int>(pmc["nskip"]) * Rcpp::as<int>(pmc["nsave"]);
  num_burnin = Rcpp::as<int>(pmc["nburn"]);
  num_thin = Rcpp::as<int>(pmc["nskip"]);
  num_display = Rcpp::as<int>(pmc["ndisplay"]);

  length_chain = num_iter/num_thin;

  saveT.set_size(length_chain, n);
  saveZ.set_size(length_chain, n);
  saveW.set_size(J, K, length_chain);
  saveXi.set_size(J, p*K, length_chain);
  saveXi0.set_size(p, K, length_chain);
  saveE.set_size(p, K*p, length_chain);
  savePsi.set_size(p, K, length_chain);
  saveG.set_size(p, K*p, length_chain);
  saveS.set_size(length_chain, K);
  saveVarphi.set_size(length_chain);
  saveOmega.set_size(p, K*p, length_chain);
  saveAlpha.set_size(p, K, length_chain);
  saveA0.set_size(length_chain);
  saveLog_py.set_size(length_chain, K);
  saveNResampled.set_size(length_chain, K);
  savePerplexity.set_size(length_chain, K);

  T = Rcpp::as<uvec>(state["t"]);

  tau_a = Rcpp::as<vec>(prior["tau_a"]);
  e0 = Rcpp::as<double>(prior["e0"]);
  E0 = Rcpp::as<mat>(prior["E0"]);
  invE0 = inv_sympd(E0);
  // iota = Rcpp::as<double>(prior["iota"]);
  a_varphi = Rcpp::as<double>(prior["a_varphi"]);
  b_varphi = Rcpp::as<double>(prior["b_varphi"]);
  merge_step = Rcpp::as<bool>(prior["merge_step"]);
  merge_par = Rcpp::as<double>(prior["merge_par"]);

  //   saveParticles = control$saveParticles
  //   if(saveParticles) {
  //     outFolder = control$outFolder
  //     if(!(outFolder %in% dir())) {
  //       dir.create(outFolder, recursive=T)
  //     }
  //     if(!('Iterations' %in% dir(outFolder))){
  //       dir.create(paste0(outFolder, '/Iterations'), recursive=T)
  //     }
  //   }
  //   verbose = control$verbose

  main_loop(state, prior);
}

void PMC::main_loop(Rcpp::List state, Rcpp::List prior) {

  // cout << "Hello world" <<endl;
  int km = 0;

  umat N(J, K);
  mat logW(J, K);
  logW.fill(log(1.0/K));
  cube xi(J, p, K);
  mat xi0(p, K);
  mat psi(p, K);
  cube G(p,p,K);
  cube E(p,p,K);
  cube Omega(p,p,K);
  mat alpha(p,K);
  rowvec z(n);
  umat S(num_particles, K);

  double a = 1;
  double a0 = 1;
  // double alpha_1 = 1;
  double a_old = a;
  double a0_old = a0;
  // double alpha_old_1 = alpha_1;
  double a_par  = sqrt(K);
  double a0_par  = sqrt(K);
  // double alpha_par_1  = sqrt(K);
  double a_count = 0;
  int a_tot = 100;
  // double alpha_count_0 = 0;
  // int alpha_tot_0 = 100;


  mat log_dQ(num_particles, 8);
  log_dQ.fill(0);
  vec log_py(K);
  vec perplexity(K);
  uvec nResampled(K);

  Rcpp::List all_particles;
  vec varphi(num_particles);
  // if (init==true) {
    cout << "initializing all particles..." << endl;
    varphi = initialVarphiParticles( prior );
    all_particles = Rcpp::as<Rcpp::List>(initialParticles( T, varphi) );
    cout << "Done" << endl;
  // }
  // if (init==false) {
  //   cout << "using fixed initial particles" << endl;
  //   varphi.fill(2.0/3.0);
  //   all_particles = Rcpp::as<Rcpp::List>(initParticles);
  // }



  for (int it=0; it<(num_iter+num_burnin); it++) {
    if((it+1)%num_display==0)
      cout << "Iteration: " << it + 1 << " of " << num_iter + num_burnin << endl;

    N.fill(0);
    for (int i=0; i<n; i++) {
      N(C(i), T(i))++;
    }

    // if(merge_step && it>0)  {
    if(merge_step && it>0)  {
      for( int k=0; k < K-1 ; k++ )  {
        if( sum(N.col(k)) > 0 )  {
          for(int kk=k+1; kk < K ; kk++)  {
            if( sum(N.col(kk)) > 0  )   {
              double kl_div = KL( xi0.col(k),
                                  xi0.col(kk),
                                  Omega.slice(k),
                                  Omega.slice(kk),
                                  alpha.col(k),
                                  alpha.col(kk) );
              if( kl_div < R::qchisq(merge_par, (double)p, 1, 0) ) {
                // cout << "k=" << k+1 << endl;
                // cout << "kk=" << kk+1 << endl;
                // cout << "xi0_k=c(" << xi0(0,k) << "," << xi0(1,k) << ");"<< endl;
                // cout << "xi0_kk=c(" << xi0(0,kk) << "," << xi0(1,kk) << ");"<< endl;
                // mat omk = Omega.slice(k);
                // cout << "Omega_k=matrix(c("<<omk(0,0)<<","<<omk(0,1)<<","<<omk(1,0)<<","<<omk(1,1)<<"),2,2)"  << endl;
                // omk = Omega.slice(kk);
                // cout << "Omega_kk=matrix(c("<<omk(0,0)<<","<<omk(0,1)<<","<<omk(1,0)<<","<<omk(1,1)<<"),2,2)"  << endl;
                // cout << "alpha_k=c(" << alpha(0,k) << "," << alpha(1,k) << ");"<< endl;
                // cout << "alpha_kk=c(" << alpha(0,kk) << "," << alpha(1,kk) << ");"<< endl;
                // cout << "KL=" << kl_div << endl;
                // cout << "Nk=" << N.col(k).t() << endl;
                // cout << "Nkk=" << N.col(kk).t() << endl;
                // cout << "N before" << N << endl;

                // Consider weighing:
                // double kw = accu(N.col(k))/(accu(N.col(k))+accu(N.col(kk)));
                // double kkw = accu(N.col(kk))/(accu(N.col(k))+accu(N.col(kk)));
                // xi0.col(k) = kw*xi0.col(k) + kkw*xi0.col(kk);

                N.col(k) = N.col(k) + N.col(kk);
                N.col(kk) = zeros<uvec>(J);
                // cout << "N after" << N << endl;
                T(find(T==kk)).fill(k);
                // cout << "merged" << endl;
                if( (it+1 > num_burnin) && ((it+1) % num_thin == 0)) {
                  cout << "Merged clusters after burn-in period (iteration " << it+1 << "). Consider longer burn-in." << endl;
                } else {
                  cout << "Merged clusters (iteration " << it+1 << ")" << endl;
                }


                // List tempSMuSigma = PriorSMuSigma(   varphi,
                //                                      Sigma_1,
                //                                      Omega_1,
                //                                      k_0,
                //                                      epsilon,
                //                                      m_1 );
                // S(kk) = Rcpp::as<unsigned>(tempSMuSigma["S"]);
                // mu.slice(kk) = Rcpp::as<mat>(tempSMuSigma["mu"]);
                // mu_0.col(kk) = Rcpp::as<vec>(tempSMuSigma["mu_0"]);
                // Omega.slice(kk) = Rcpp::as<mat>(tempSMuSigma["Omega"]);
                // Sigma.slice(kk) = Rcpp::as<mat>(tempSMuSigma["Sigma"]);
              }
            }
          }
        }
      }
    }


    sampleA0(a0, N, a_par);
    // cout << "did A0" << endl;

    a0_old = a0;
    a0 = sampleA0( a0, N, a_par );
    if( it <= num_burnin ) {
      if( a0 != a0_old )
        a_count++;

      if( (it+1)  % a_tot == 0) {
        if( a_count < 30 )
          a_par *= 1.1;
        if( a_count > 50 )
          a_par *= 0.9;
        a_count = 0;
      }
    } else {
      if( a0 != a0_old )
        a_count++;
    }

    logW = sampleLogWs(N, a0);
    // cout << "did W" << endl;
    Rcpp:List varphiList = sampleVarphi( S, prior );
    varphi = Rcpp::as<vec>(varphiList["varphi"]);
    log_dQ.col(0) = Rcpp::as<vec>(varphiList["log_dq"]);
    // cout << "did varphi" << endl;

    for (int k=0; k<K; k++) {
      // cout << "k= " << k+1 << endl;
      // if (accu(N.col(k))==0) cout << "start empty iter" << endl;
      Rcpp::List iterSummary = iter( T, k, N, all_particles[k], log_dQ, varphi, prior);
      // if (accu(N.col(k))==0) cout << "end empty iter" << endl;

      all_particles[k] = Rcpp::as<Rcpp::List>(iterSummary["particles"]);

      Rcpp::List temp = all_particles[k];

      xi.slice(k) = mean(Rcpp::as<cube>(temp["xi"]), 2);
      xi0.col(k) = mean(Rcpp::as<mat>(temp["xi0"]),0).t();
      psi.col(k) = mean(Rcpp::as<mat>(temp["psi"]),0).t();
      G.slice(k) = reshape(mean(Rcpp::as<mat>(temp["G"]),0), p ,p);
      E.slice(k) = reshape(mean(Rcpp::as<mat>(temp["E"]),0), p ,p);
      z.cols(find(T==k)) = mean(Rcpp::as<mat>(temp["z"]),0);
      S.col(k) = Rcpp::as<uvec>(temp["S"]);

      // if (accu(N.col(k)==0)) cout << "empty trans params" << endl;
      Omega.slice(k) = G.slice(k) + psi.col(k) * psi.col(k).t();
      mat inv_Omega = inv_sympd(Omega.slice(k));
      vec numerator = arma::sqrt(diagmat(Omega.slice(k).diag())) * inv_Omega * psi.col(k);
      double denominator = as_scalar(arma::sqrt(1 - psi.col(k).t() * inv_Omega * psi.col(k)) );
      alpha.col(k) = numerator/denominator;
      // if (accu(N.col(k)==0))  cout << "empty end trans params" << endl;

      log_py(k) = Rcpp::as<double>(iterSummary["log_py"]);
      perplexity(k) = Rcpp::as<double>(iterSummary["perplexity"]);
      nResampled(k) = Rcpp::as<int>(iterSummary["nResampled"]);
    }

    T = sampleT( xi, Omega, alpha, logW );

    if( (it+1 > num_burnin) && ((it+1) % num_thin == 0))
    {
      saveT.row(km) = T.t();
      saveZ.row(km) = z;
      saveW.slice(km) = exp(logW);
      saveXi.slice(km) = reshape( mat(xi.memptr(), xi.n_elem, 1, false), J, K*p);
      saveXi0.slice(km) = xi0;
      savePsi.slice(km) = psi;
      saveG.slice(km) = reshape( mat(G.memptr(), G.n_elem, 1, false), p, K*p);
      saveE.slice(km) = reshape( mat(E.memptr(), E.n_elem, 1, false), p, K*p);
      saveLog_py.row(km) = log_py.t();
      saveS.row(km) = (sum(S,0)>K/2);
      saveVarphi(km) = mean( varphi );
      saveA0(km) = a0;
      savePerplexity.row(km) = perplexity.t();
      saveNResampled.row(km) = nResampled.t();
      // saveSigma.slice(km) = reshape( mat(Sigma.memptr(), Sigma.n_elem, 1, false), p, K*p);
      saveOmega.slice(km) = reshape( mat(Omega.memptr(), Omega.n_elem, 1, false), p, K*p);
      saveAlpha.slice(km) = alpha;
      km++;
    }
  //
  } // end it loop
}

//--------------------------------------------------------

arma::mat PMC::sampleLogWs(   arma::umat N,
                               double a0) {
  mat logW(J,K);
  for(int j=0; j<J; j++) {
    logW.row(j) = rDirichlet( N.row(j).t() +  a0 * ones<vec>(K) / K ).t();
  }

  return logW;
}


Rcpp::List PMC::initialParticles( arma::uvec T, arma::vec varphi ) {
  // Given the arguments, initialPoints creates a population
  // of particles from simple proposals: sample z iid N(0,1).
  // If N_k>0:
  // psi_k as mean of the normal part of the full conditional (with xi0k=mean(Y_k)),
  // xi_jk, xi0k as mean of full conditional [xi0k|psi_k, Sk=0, ...],
  // G as mean of IW part of full conditional (taking Lambda=0pxp, Sk)
  // If N_k==0:
  // psi as mean of the normal part of the full conditional
  // (with all data and xi0k=mean(Y)),
  // xi0k as mean of full conditional with full data, latter psi
  // half of xi_j as above xi0k, half as xi0k + noise
  // G as mean of IW part of full conditional (taking Lambda=0pxp, S)
  //
  // E_k from E_k|S_k=0,...  which is the prior for E_k
  // S_k and varphi for now from their priors
  // *** consider sampling S_k from prior and then the rest accordignly (instead of all from S_k=0)

  Rcpp::List all_particles(K);

  for ( int k=0; k<K; k++ ) {
    uvec T_k = arma::find(T==k);
    int n_k = T_k.n_elem;

    mat z_k;
    mat psi_k(num_particles,p);
    psi_k.fill(0);
    cube xi_k(J, p, num_particles);
    mat xi_0k(num_particles, p);
    mat G_k(num_particles, pow(p,2));
    mat E_k(num_particles, pow(p,2));
    uvec S_k(num_particles);

    if( n_k > 0 ) {
      z_k = randn<mat>( num_particles, n_k );
      mat Y_k = Y.rows(T_k);
      rowvec mean_y = mean(Y_k, 0);
      mat absz = abs(z_k);
      vec mean_absz = mean(absz, 1);

      mat zmat = absz.each_col() - mean_absz;
      mat ymat = Y_k.each_row() - mean_y;
      vec devz = sum(pow(zmat,2), 1);
      for (int iN=0; iN<num_particles; iN++) {
        for (int ip=0; ip<p; ip++) {
          psi_k(iN, ip) = as_scalar(accu(zmat.row(iN).t() % ymat.col(ip)) / devz(iN));
        }
      }

      xi_0k = -1.0*(psi_k % repmat(mean_absz, 1, p));
      xi_0k.each_row() += mean_y;

      uvec C_k = C(T_k);
      for (int j=0; j<J; j++) {
        uvec j_indices = find(C_k==j);
        int n_jk = j_indices.n_elem;
        if (n_jk>0) {
          mat z_jk = z_k.cols(j_indices);
          mat Y_jk = Y_k.rows(j_indices);
          rowvec mean_yjk = mean(Y_jk, 0);
          mat absz_jk = abs(z_jk);
          vec mean_absz_jk = mean(absz_jk, 1);

          mat xi_jk = -1.0*(psi_k % repmat(mean_absz_jk, 1, p));
          xi_jk.each_row() += mean_yjk;
          xi_k.subcube(j,0,0,j,p-1,num_particles-1) = xi_jk.t();
        } else {
          xi_k.subcube(j,0,0,j,p-1,num_particles-1) = xi_0k.t();
        }
      }

      for ( int iN=0; iN<num_particles; iN++ ) {
        mat e;
        int nG;
        if ( n_k > (uint)p + 2 ) {
          e = ( Y_k.each_row() - xi_0k.row(iN) );
          nG = n_k;
        } else {
          e = ( Y.each_row() - xi_0k.row(iN) );
          nG = T.n_elem;
        }
        e -= ( repmat(psi_k.row(iN), nG, 1) % repmat(abs(z_k.row(iN).t()), 1, p) );
        mat Gsum(p, p);
        Gsum.fill(0);
        for (int i=0; i<nG; i++ ) {
          Gsum += e.row(i).t() * e.row(i);
        }
        G_k.row(iN) = vectorise(Gsum / nG).t();
      }
    }
    else { // if N_k==0, do same with all data, but return empty z
      int n_all = Y.n_rows;

      mat temp_z_k = randn<mat>( num_particles, n_all );
      rowvec mean_y = mean(Y, 0);
      mat absz = abs(temp_z_k);
      vec mean_absz = mean(absz, 1);

      mat zmat = absz.each_col() - mean_absz;
      mat ymat = Y.each_row() - mean_y;
      vec devz = sum(pow(zmat,2), 1);

      for (int iN=0; iN<num_particles; iN++) {
        for (int ip=0; ip<p; ip++) {
          psi_k(iN, ip) = as_scalar(accu(zmat.row(iN).t() % ymat.col(ip)) / devz(iN));
        }
      }

      xi_0k = -1.0*(psi_k % repmat(mean_absz, 1, p));
      xi_0k.each_row() += mean_y;

      vec u = arma::randu<vec>(J);
      double priorVarphi = a_varphi/(a_varphi+b_varphi);
      for (int j=0; j<J; j++) {
        if ( u(j) < priorVarphi ) {
          xi_k.subcube(j,0,0,j,p-1,num_particles-1) = trans( xi_0k + mvrnormArma(num_particles, mean(xi_0k, 0).t(), E0/(e0 - p - 1.0)) );
        } else {
          xi_k.subcube(j,0,0,j,p-1,num_particles-1) = xi_0k.t();
        }
      }

      for ( int iN=0; iN<num_particles; iN++ ) {
        mat e = ( Y.each_row() - xi_0k.row(iN) );
        e -= ( repmat(psi_k.row(iN), n_all, 1) % repmat(abs(temp_z_k.row(iN).t()), 1, p) );
        mat Gsum(p, p);
        Gsum.fill(0);
        for (int i=0; i<n_all; i++ ) {
          Gsum += e.row(i).t() * e.row(i);
        }
        G_k.row(iN) = vectorise(Gsum / n_all).t();
      }
      z_k.set_size( num_particles, 0 );
    } // end N_k==0

    S_k.fill(0);
    for ( int iN=0; iN<num_particles; iN++ ) {
      mat E = inv_sympd(rWishartArma(inv_sympd(E0), e0 ));
      E_k.row(iN) = vectorise(E).t();
      // varphi(iN) = R::rbeta(a_varphi, b_varphi);
      double u = R::runif(0,1);
      if (u<varphi(iN)) S_k(iN) = 1;
    }

    Rcpp::List particles = Rcpp::List::create(
      Rcpp::Named( "z" ) = z_k,
      Rcpp::Named( "psi" ) = psi_k,
      Rcpp::Named( "xi" ) = xi_k,
      Rcpp::Named( "xi0" ) = xi_0k,
      Rcpp::Named( "G" ) = G_k,
      Rcpp::Named( "E" ) = E_k,
      Rcpp::Named( "S" ) = S_k);

    all_particles(k) = particles;
  }
  return all_particles;
}

vec PMC::initialVarphiParticles( Rcpp::List prior ) {

  vec varphi(num_particles);
  double a_varphi = Rcpp::as<double>(prior["a_varphi"]);
  double b_varphi = Rcpp::as<double>(prior["b_varphi"]);

  for ( int iN=0; iN<num_particles; iN++ )
    varphi(iN) = R::rbeta(a_varphi, b_varphi);

  return varphi;
}


Rcpp::List PMC::sampleXi(mat Y_k, uvec C_k, uvec N_k, Rcpp::List particles ) {
  // Given the arguments, this function returns a population of MC draws for the values of the variable Xi, in the p-variate skew-N model.
  int nk = Y_k.n_rows;

  uvec S = Rcpp::as<umat>(particles["S"]);
  mat z = Rcpp::as<mat>(particles["z"]);
  mat psi = Rcpp::as<mat>(particles["psi"]);
  mat G = Rcpp::as<mat>(particles["G"]);
  mat E = Rcpp::as<mat>(particles["E"]);
  mat xi0 = Rcpp::as<mat>(particles["xi0"]);

  cube xi(J, p, num_particles);
  xi.fill(0);
  mat log_dxi(num_particles, J);

  if(nk == 0) {
    // cout << "empty xi" << endl;
    for ( int iN=0; iN<num_particles; iN++ ) {
      if ( S(iN)==0 ) {
        xi.slice(iN) = repmat( xi0.row(iN), J, 1 );
        log_dxi.row(iN).fill(0);
      }
      else { // if S(iN)==1
        mat tempE = reshape(E.row(iN), p, p);
        double sgn, ld;
        log_det( ld, sgn, tempE );
        vec eigval; mat eigvec;
        eig_sym(eigval , eigvec,  tempE) ;
        mat invSqrtL_invU = diagmat(1.0/sqrt(eigval)) * inv(eigvec);
        for(int j=0; j<J; j++) {
          xi.slice(iN).row(j) = trans( mvnrnd(xi0.row(iN).t(), tempE) );
          vec normedXi = invSqrtL_invU * trans(xi.slice(iN).row(j) - xi0.row(iN));
          log_dxi(iN, j) = as_scalar(-log2pi * (p/2.0) -0.5*normedXi.t()*normedXi -0.5*sgn*ld );
          // log_dxi(iN, j) = as_scalar(dmvnrm_arma_precision(
          //                  xi.slice(iN).row(j), xi0.row(iN),
          //                  inv_sympd(tempE) ) );
        } // end j loop
      } // end S(iN==1)|nk=0
    } // end iN loop
    // cout << "end empty xi" << endl;
  } else { // if nk>0:
    for ( int iN=0; iN<num_particles; iN++ ) {
      if ( S(iN)==0 ) {
        xi.slice(iN) = repmat( xi0.row(iN), J, 1 );
        log_dxi.row(iN).fill(0);
      } else { // if S(iN)==1
        mat invE = inv_sympd(reshape(E.row(iN),p,p));
        mat invG = inv_sympd(reshape(G.row(iN),p,p));
        for(int j=0; j<J; j++) {
          int njk = N_k(j);
          if ( njk ==0 ) {
            // xi.slice(iN).row(j) = xi0.row(iN);
            xi.slice(iN).row(j) = trans( mvnrnd(xi0.row(iN).t(), reshape(E.row(iN),p,p)) );
            // log_dxi(iN, j) = 0;
            log_dxi(iN, j) = as_scalar(dmvnrm_arma_precision(
                                              xi.slice(iN).row(j),
                                              xi0.row(iN),
                                              reshape(E.row(iN),p,p)) );
          } else { // if njk>0
            uvec jk_idx = find(C_k==j);
            mat Yjk = Y_k.rows(jk_idx);
            mat zin = z.row(iN);
            rowvec zjk = zin.cols(jk_idx);
            mat v = inv_sympd(invE + njk*invG);
            rowvec mean_y = mean(Yjk, 0);
            vec m = invE * xi0.row(iN).t() + njk*invG * trans((mean_y - (psi.row(iN) * mean(abs(zjk)))));
            xi.slice(iN).row(j) = trans(mvnrnd(v * m, v));
            log_dxi(iN, j) = as_scalar(dmvnrm_arma_precision(
                                  xi.slice(iN).row(j), trans(v * m),
                                  invE + njk*invG ) );
          }
        } // end j loop
      } // end S(iN)==1 if
    } // end iN (particles) loop
  } // end nk>0 loop

  vec vec_log_dxi = sum(log_dxi, 1);

  return Rcpp::List::create(
    Rcpp::Named( "xi" ) = xi,
    Rcpp::Named( "log_dq" ) = vec_log_dxi
  );
}

Rcpp::List PMC::sampleG( mat Y_k, uvec C_k, uvec N_k,
                         Rcpp::List particles, Rcpp::List prior ) {

  int nk = Y_k.n_rows;

  uvec S = Rcpp::as<umat>(particles["S"]);
  cube xi = Rcpp::as<cube>(particles["xi"]);
  mat xi0 = Rcpp::as<mat>(particles["xi0"]);
  mat absz = abs(Rcpp::as<mat>(particles["z"]));
  mat psi = Rcpp::as<mat>(particles["psi"]);
  double m = Rcpp::as<double>(prior["m0"]);
  mat Lambda = Rcpp::as<mat>(prior["Lambda"]);

  mat G(num_particles, pow(p,2));
  vec log_dG(num_particles);
  int iN;
  mat ek(nk,p);
  mat zp(nk,p);
  mat Lambda_k(p,p);
  mat g(p,p);

  if (nk==0) {
    // cout << "empty G" << endl;
    for ( iN=0; iN<num_particles; iN++) {
      g = inv_sympd(rWishartArma(inv_sympd(Lambda), m ));
      G.row(iN) = vectorise(g).t();
      log_dG(iN) = dIWishartArma(g, m, Lambda);
    }
    // cout << "end empty G" << endl;
  } else { // if nk>0
    for ( iN=0; iN<num_particles; iN++ ) {
      mat Sk(p,p);
      Sk.fill(0);
      if ( S(iN)==0 ) {
        // cout << "at S==0" << endl;
        // #pragma omp parallel for private(ip, e, zp, C_k, g)
        ek = Y_k.each_row() - xi0.row(iN);
        zp = repmat(absz.row(iN).t(), 1 , p);
        zp.each_row() %= psi.row(iN);
        ek -= zp;
        Sk = ek.t() * ek;
      } else { // if S(iN)==1
        // cout << "at S==1" << endl;
        for (int j=0; j<J; j++) {
          int njk = N_k(j);
          if ( njk > 0 ) {
            uvec jk_idx = find(C_k==j);
            mat ejk(njk, p);
            mat zpj(njk, p);
            rowvec xi_j = xi.subcube(j, 0, iN, j, p-1, iN);
            mat Y_jk = Y_k.rows(jk_idx);
            ejk =  Y_jk.each_row() - xi_j;
            vec zrow = trans( absz.row(iN) );
            zpj = repmat( zrow(jk_idx) , 1 ,p);
            zpj.each_row() %= psi.row(iN);
            ejk -= zpj;
            Sk += ejk.t() * ejk;
          }
        }
      }
      if ( accu(abs(Lambda))==0 )
        Lambda_k = Sk;
      else
        Lambda_k = inv_sympd(Lambda) + Sk;
      g = inv_sympd(rWishartArma(inv_sympd(Lambda_k), nk+m ));
      G.row(iN) = vectorise(g).t();
      log_dG(iN) = dIWishartArma(g, nk+m, Lambda_k);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named( "G" ) = G,
    Rcpp::Named( "log_dq" ) = log_dG  );
}


Rcpp::List PMC::samplePsi( mat Y_k, uvec C_k, uvec N_k, Rcpp::List particles,
                           Rcpp::List prior) {

  int nk = Y_k.n_rows;
  int p = Y_k.n_cols;

  uvec S = Rcpp::as<umat>(particles["S"]);
  cube xi = Rcpp::as<cube>(particles["xi"]);
  mat G = Rcpp::as<mat>(particles["G"]);
  mat xi0 = Rcpp::as<mat>(particles["xi0"]);
  mat absz = abs(Rcpp::as<mat>(particles["z"]));
  vec sums_z2 = sum(pow(absz,2), 1);

  mat psi(num_particles, p);
  vec log_dpsi(num_particles);

  if(nk == 0) {
    // cout << "empty psi" << endl;
    // add sampling from prior on n-ball
    for ( int iN=0; iN<num_particles; iN++ ) {
      vec x = randn(p);
      double r = norm(x);
      double u = R::runif(0,1);
      // sample from Sigma prior (take psi*psi element here to be zero)
      double m = Rcpp::as<double>(prior["m0"]);
      mat Lambda = Rcpp::as<mat>(prior["Lambda"]);
      mat g = inv_sympd(rWishartArma(inv_sympd(Lambda), m ));
      rowvec h = sqrt(g.diag().t());
      mat invD = inv(diagmat(h));
      double detOmega =  det(invD * g * invD);

      vec delta = pow(u, 1/p) * (1 / r) * pow(detOmega, p/2) * x;
      psi.row(iN) = trans( diagmat(h) * delta );
      double logSphere = (p/2.0) * log(M_PI) - lgamma(p/2.0+1.0);
      log_dpsi(iN) = 1 / (logSphere * log(detOmega));
    }
    // cout << "end empty psi" << endl;
  } else { // if nk>0:
    for ( int iN=0; iN<num_particles; iN++ ) {
      rowvec mpsi(p);
      mpsi.fill(0);
      if ( S(iN)==0 ) {
        mat Y_dmean = Y_k.each_row() - xi0.row(iN);
        Y_dmean.each_col() %= absz.row(iN).t();
        mpsi = sum(Y_dmean, 0) / sums_z2(iN);
      } else { // if S(iN)==1
        for (int j=0; j<J; j++) {
          int njk = N_k(j);
          if ( njk > 0 ) {
            uvec jk_idx = find(C_k==j);
            rowvec xi_j = xi.subcube(j, 0, iN, j, p-1, iN);
            mat Y_jk = Y_k.rows(jk_idx);
            mat Yjk_dmean = Y_jk.each_row() - xi_j;
            vec zrow = trans( absz.row(iN) );
            Yjk_dmean.each_col() %= zrow(jk_idx);
            mpsi += sum(Yjk_dmean, 0);
          }
        } // end j loop
        mpsi /= sums_z2(iN);
      }
      mat vpsi = reshape(G.row(iN) / sums_z2(iN), p, p);
      psi.row(iN) = (arma::mvnrnd(mpsi.t(), vpsi)).t();
      log_dpsi(iN) = as_scalar(
        dmvnrm_arma_precision(  psi.row(iN),
                                mpsi,
                                inv_sympd(vpsi) ));
    } // end iN (particles) loop
  } // end nk>0 if

  return Rcpp::List::create(
    Rcpp::Named( "psi" ) = psi,
    Rcpp::Named( "log_dq" ) = log_dpsi
  );

}

Rcpp::List PMC::sampleZ( mat Y_k, uvec C_k, Rcpp::List particles ) {

  int nk = Y_k.n_rows;
  vec log_dZ(num_particles);
  mat z(num_particles, nk);

  if ( nk==0 ) {
    // cout << "empty z" << endl;
    log_dZ.fill(0);
    // cout << "end empty z" << endl;
  } else {
    uvec S = Rcpp::as<umat>(particles["S"]);
    cube xi = Rcpp::as<cube>(particles["xi"]);
    mat xi0 = Rcpp::as<mat>(particles["xi0"]);
    mat psi = Rcpp::as<mat>(particles["psi"]);
    mat G = Rcpp::as<mat>(particles["G"]);
    int i;
    double m;
    mat xi_j;
    double sgn, absz;
    vec tp(nk); vec ldZ(nk); vec u(nk);

    for (int ip=0; ip<num_particles; ip++) {
      mat invG = inv_sympd(reshape(G.row(ip), p, p));
      double v = 1.0/(1.0 + as_scalar(psi.row(ip) * invG * psi.row(ip).t()));
      double sv = sqrt(v);
      u = randu(nk);
      for ( i=0; i<nk; i++ ) {
        if (S(ip)==0) {
          m = as_scalar(v * (psi.row(ip) * invG * (Y_k.row(i) - xi0.row(ip)).t()));
        } else {
          xi_j = xi.subcube(C_k(i), 0, ip, C_k(i), p-1, ip);
          m = as_scalar(v * (psi.row(ip) * invG * (Y_k.row(i) - xi_j).t()));
        }
        if (u(i) < 0.5) {
          sgn = -1.0;
        } else {
          sgn = 1.0;
        }
        absz = rtruncnormArma(m, sv, 0);
        z(ip, i) = sgn * absz;
        tp(i) = R::pnorm(0, -m, sv, true, true);
        ldZ(i) = R::dnorm(absz, m, sv, true);
      }
      log_dZ(ip) = accu( ldZ - tp - log(2) );
    }
  }

  return Rcpp::List::create(
    Rcpp::Named( "z" ) = z,
    Rcpp::Named( "log_dq" ) = log_dZ
  );
}


Rcpp::List PMC::sampleXi0(mat Y_k, uvec N_k, Rcpp::List particles ) {
  int nk = Y_k.n_rows;

  uvec S = Rcpp::as<umat>(particles["S"]);
  mat psi = Rcpp::as<mat>(particles["psi"]);
  mat G = Rcpp::as<mat>(particles["G"]);
  mat E = Rcpp::as<mat>(particles["E"]);
  cube xi = Rcpp::as<cube>(particles["xi"]);
  mat absz = abs(Rcpp::as<mat>(particles["z"]));

  mat xi0(num_particles, p);
  xi0.fill(0);
  vec log_dxi0(num_particles);

  if(nk == 0) {
    // sample from ...
    // cout << "empty xi0" << endl;
    rowvec mxi0 = mean(Y,0);
    mat vxi0 = 100*cov(Y);
    xi0 = mvrnormArma(num_particles, mxi0.t(), vxi0);
    log_dxi0 =  dmvnrm_arma_precision(xi0,
                                      mxi0,
                                      inv_sympd(vxi0) );
    // cout << "end empty xi0" << endl;
  } else { // if nk>0:
    for ( int iN=0; iN<num_particles; iN++ ) {
      rowvec mxi0(p);
      mat vxi0(p,p);
      if ( S(iN)==0 ) {
        rowvec mean_y = mean(Y_k, 0);
        double mean_absz = mean(absz.row(iN));
        mxi0 = mean_y - psi.row(iN) * mean_absz;
        vxi0 = reshape(G.row(iN), p, p)/nk;
      } else { // if S(iN)==1
        uvec jk_idx = find(N_k>0);
        mat xi_slice = xi.slice(iN);
        mxi0 = mean(xi_slice.rows(jk_idx), 0);
        vxi0 = reshape(E.row(iN), p, p)/jk_idx.n_elem;
      } // end S(iN)==1 if
      xi0.row(iN) = mvrnormArma(1, mxi0.t(), vxi0);
      log_dxi0(iN) =  as_scalar(dmvnrm_arma_precision(
                                  xi0.row(iN),
                                  mxi0,
                                  inv_sympd(vxi0) ) );
    } // end iN (particles) loop
  } // end nk>0 loop

  return Rcpp::List::create(
    Rcpp::Named( "xi0" ) = xi0,
    Rcpp::Named( "log_dq" ) = log_dxi0
  );
}

Rcpp::List PMC::sampleE( uvec N_k, Rcpp::List particles, Rcpp::List prior ) {

  // cout << "e0=" << e0 << endl;
  int nk = accu(N_k);

  uvec S = Rcpp::as<umat>(particles["S"]);
  cube xi = Rcpp::as<cube>(particles["xi"]);
  mat xi0 = Rcpp::as<mat>(particles["xi0"]);
  // double e0 = Rcpp::as<double>(prior["e0"]);
  // mat E0 = Rcpp::as<mat>(prior["E0"]);

  mat E(num_particles, pow(p,2));
  vec log_dE(num_particles);
  int iN;

  if (nk==0) {
    // cout << "empty E" << endl;
    for ( iN=0; iN<num_particles; iN++) {
      // mat e = inv_sympd(rWishartArma( inv_sympd(E0), e0 ));
      mat e = inv_sympd(rWishartArma( invE0, e0 ));
      E.row(iN) = vectorise(e).t();
      log_dE(iN) = dIWishartArma(e, e0, E0);
    }
    // cout << "end empty E" << endl;
  } else { // if nk>0
    for ( iN=0; iN<num_particles; iN++ ) {
      mat Sk(p,p);
      Sk.fill(0);
      if ( S(iN)==0 ) {
        mat e = inv_sympd(rWishartArma( invE0, e0 ));
        E.row(iN) = vectorise(e).t();
        log_dE(iN) = dIWishartArma(e, e0, E0);
      } else { // if S(iN)==1
        uvec jk_idx = find(N_k>0);
        double ek = e0 + jk_idx.n_elem;
        mat xi_slice = xi.slice(iN).rows(jk_idx);
        mat xi_dmean = xi_slice.each_row() - xi0.row(iN);
        mat mE = ( invE0 + xi_dmean.t() * xi_dmean );
        mat e = inv( rWishartArma( mE, ek ) );
        E.row(iN) = vectorise(e).t();
        log_dE(iN) = dIWishartArma(e, ek, inv_sympd(mE));
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named( "E" ) = E,
    Rcpp::Named( "log_dq" ) = log_dE  );
}

Rcpp::List PMC::sampleS( mat Y_k, uvec C_k, uvec N_k,
                         Rcpp::List particles, vec varphi) {

  int nk = accu(N_k);

  cube xi = Rcpp::as<cube>(particles["xi"]);
  mat xi0 = Rcpp::as<mat>(particles["xi0"]);
  mat psi = Rcpp::as<mat>(particles["psi"]);
  mat G = Rcpp::as<mat>(particles["G"]);
  mat E = Rcpp::as<mat>(particles["E"]);
  mat z = Rcpp::as<mat>(particles["z"]);
  // vec varphi = Rcpp::as<vec>(particles["varphi"]);

  uvec S(num_particles);
  S.fill(0);
  vec log_dS(num_particles);
  int iN;

  if (nk==0) { // Sample from prior
    // cout << "empty S" << endl;
    for ( iN=0; iN<num_particles; iN++) {
      double u = R::runif(0,1);
      if (u<varphi(iN)) {
        S(iN) = 1;
        log_dS(iN) = log(varphi(iN));
      } else
        log_dS(iN) = log(1-varphi(iN));
    }
    // cout << "end empty S" << endl;
  } else { // if nk>0
    for ( iN=0; iN<num_particles; iN++ ) {
      mat g = reshape(G.row(iN), p, p);
      double lnmr = 0;
      double ldnm = 0;
      // likelihood
      mat mxi;
      for (int i=0; i<nk; i++) {
        rowvec xi_S1 = xi.subcube(C_k(i), 0, iN, C_k(i), p-1, iN);
        lnmr += as_scalar( dmvnrm_arma_precision( Y_k.row(i),
                                                  xi_S1 + psi.row(iN)*abs(z(i)),
                                                  g ) );
        rowvec xi_S0 = xi0.row(iN);
        ldnm += as_scalar( dmvnrm_arma_precision( Y_k.row(i),
                                                  xi_S0 + psi.row(iN)*abs(z(i)),
                                                  g ) );
      }
      // add prior on xi_jk
      mat xi_jk = xi.subcube(0, 0, iN, J-1, p-1, iN);
      lnmr += accu( dmvnrm_arma_precision(  xi_jk,
                                            xi0.row(iN),
                                            reshape( E.row(iN), p, p ) )  );
      lnmr += log(varphi(iN));
      ldnm += log(1-varphi(iN));
      double probSk1 = exp( lnmr - log_exp_x_plus_exp_y( lnmr, ldnm )  );
      double u = R::runif(0,1);
      if ( u < probSk1 ) {
        S(iN) = 1;
        log_dS(iN) = log(probSk1);
      } else
        log_dS(iN) = log(1-probSk1);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named( "S" ) = S,
    Rcpp::Named( "log_dq" ) = log_dS  );
}


Rcpp::List PMC::sampleVarphi( umat S, Rcpp::List prior ) {
  vec varphi(num_particles);
  vec log_dvphi(num_particles);
  double a_varphi = Rcpp::as<double>(prior["a_varphi"]);
  double b_varphi = Rcpp::as<double>(prior["b_varphi"]);

  for (int iN=0; iN<num_particles; iN++) {
    uvec nzero = find( S.row(iN)==0 );
    int s0 = nzero.n_elem;
    int s1 = K - s0;
    varphi(iN) = as<double>(rbeta(1, a_varphi + s1, b_varphi + s0));
    log_dvphi(iN) = dBeta(varphi(iN), a_varphi + s1, b_varphi + s0);
  }
  return Rcpp::List::create(
    Rcpp::Named( "varphi" ) = varphi,
    Rcpp::Named( "log_dq" ) = log_dvphi  );
}

double PMC::sampleA0(double a0, arma::umat N, double a_par)
{
  double output = a0;
  double log_ratio = 0;
  double temp = rgammaBayes(  pow( a0, 2 ) * a_par,
                              a0 * a_par );

  log_ratio += R::dgamma(a0, pow(temp,2)* a_par, 1/temp/a_par, 1);
  log_ratio -= R::dgamma(temp, pow(a0,2)* a_par, 1/a0/a_par, 1);
  log_ratio += R::dgamma(temp, tau_a(0), 1/tau_a(1), 1);
  log_ratio -= R::dgamma(a0, tau_a(0), 1/tau_a(1), 1);

  for(int j = 0; j < J; j++)  {
    log_ratio += marginalLikeDirichlet( N.row(j).t(), (temp/K)*ones<vec>(K)  );
    log_ratio -= marginalLikeDirichlet( N.row(j).t(), (a0/K)*ones<vec>(K)  );
  }
  if( exp(log_ratio) > R::runif(0,1) )
    output = temp;

  return output;
}

// ----------------------------------------------------------------------------------

arma::vec PMC::logPriorDens( Rcpp::List particles, vec varphi,
                             Rcpp::List prior ) {

  uvec S = Rcpp::as<uvec>(particles["S"]);
  cube xi = Rcpp::as<cube>(particles["xi"]);
  mat xi0 = Rcpp::as<mat>(particles["xi0"]);
  mat psi = Rcpp::as<mat>(particles["psi"]);
  mat G = Rcpp::as<mat>(particles["G"]);
  mat E = Rcpp::as<mat>(particles["E"]);

  a_varphi = Rcpp::as<double>(prior["a_varphi"]);
  b_varphi = Rcpp::as<double>(prior["b_varphi"]);
  double m = Rcpp::as<double>(prior["m0"]);
  mat Lambda = Rcpp::as<mat>(prior["Lambda"]);

  mat h(num_particles, p);
  vec logDetOmega(num_particles);
  vec logDetSigma(num_particles);
  vec logPriorG(num_particles);
  vec logPriorE(num_particles);
  vec logPriorS = log(1 - varphi);
  logPriorS( find(S==1) ) = log( varphi( find(S==1) ) );
  vec logPriorVarphi(num_particles);

  int iN;
  mat Sigma(p,p);
  mat invD(p,p);
  double logpxi = 0;

// #pragma omp parallel for private(ip, Sigma, invD)
  for ( iN=0; iN<num_particles; iN++ ) {
    Sigma = reshape(G.row(iN), p, p) + psi.row(iN).t() * psi.row(iN);
    // if(!all(eigen(Sigma.iN)$values>0)) out[iN] = 1
    logDetSigma(iN) = log(det(Sigma));
    h.row(iN) = sqrt(Sigma.diag().t());
    mat invD = inv(diagmat(h.row(iN)));
    logDetOmega(iN) =  log(det(invD * Sigma * invD));

    mat tempE = reshape(E.row(iN), p, p);
    // temppriore(iN) = log(det(tempE));
    if (S(iN)==1) {
      for (int j=0; j<J; j++) {
        logpxi += as_scalar( dmvnrm_arma_precision(
                           xi.slice(iN).row(j), xi0.row(iN),
                           inv_sympd(tempE) ));
      }
    } else {
      logpxi = 0;
    }

    logPriorE(iN) = dIWishartArma(tempE, e0, E0);
    logPriorVarphi(iN) = dBeta( varphi(iN), a_varphi, b_varphi );
    if (m==0) {
      -((p+1.0)/2.0) * logDetSigma(iN);
    } else {
      logPriorG(iN) = dIWishartArma( Sigma, m, Lambda);
    }
  }
  // pos = 1 - out

  double logSphere = (p/2.0) * log(M_PI) - lgamma(p/2.0+1.0);

  vec logPriorDens = logpxi - logSphere - 0.5 * logDetOmega +
                        logPriorG + logPriorE + logPriorS + logPriorVarphi -
                        sum(log(h), 1);

  //  logPriorDens[which(particlesStar$pos == 0)] = rep(-Inf, sum(particlesStar$pos == 0))
  return(logPriorDens);
}

arma::vec PMC::logPostDens( mat Y_k, uvec C_k, uvec N_k,
                            Rcpp::List particles, vec varphi,
                            Rcpp::List prior ) {

  int nk = Y_k.n_rows;
  vec loglikelihood(num_particles);
  loglikelihood.fill(0);

  if (nk>0) {
      uvec S = Rcpp::as<umat>(particles["S"]);
      cube xi = Rcpp::as<cube>(particles["xi"]);
      mat xi0 = Rcpp::as<mat>(particles["xi0"]);
      mat psi = Rcpp::as<mat>(particles["psi"]);
      mat G = Rcpp::as<mat>(particles["G"]);
      mat absz = abs(Rcpp::as<mat>(particles["z"]));

      vec sums_z2(num_particles);
      int iN;
      // mat e(nk,p);
      mat zp;
      mat g(p,p);
      mat ee(num_particles, pow(p,2));
      vec detG(num_particles);
      vec PsiVec(num_particles);

    // #pragma omp parallel for private(ip, e, zp, Lambda, g)
      for (iN=0; iN<num_particles; iN++) {
        if (S(iN)==0) {
          sums_z2 = arma::sum(absz%absz, 1);
          mat e = Y_k.each_row() - xi0.row(iN);
          zp = repmat(absz.row(iN).t(), 1, p);
          zp.each_row() %= psi.row(iN);
          e -= zp;
          ee = e.t() * e;
          g = reshape(G.row(iN), p, p);
          detG(iN) = det(g);
          PsiVec(iN) = accu(inv(g) % ee);
          loglikelihood(iN) = as_scalar((- nk/2.0) * log(detG(iN)) - 0.5 * PsiVec(iN) + (- 0.5) * sums_z2.row(iN));
        } else { // if S(k)==1
          for (int j=0; j<J; j++) {
            int njk = N_k(j);
            if ( njk > 0 ) {
              uvec jk_idx = find(C_k==j);
              mat absz_jk = absz.cols(jk_idx);
              sums_z2 = arma::sum(absz_jk % absz_jk, 1);
              mat Y_jk = Y_k.rows(jk_idx);
              rowvec xi_j = xi.subcube(j, 0, iN, j, p-1, iN);
              mat e = Y_jk.each_row() - xi_j;
              zp = repmat(absz_jk.row(iN).t(), 1, p);
              zp.each_row() %= psi.row(iN);
              e -= zp;
              ee = e.t() * e;
              g = reshape(G.row(iN), p, p);
              detG(iN) = det(g);
              PsiVec(iN) = accu(inv(g) % ee);
              loglikelihood(iN) += as_scalar((- njk/2.0) * log(detG(iN)) - 0.5 * PsiVec(iN) + (- 0.5) * sums_z2.row(iN));
            }
          }
        }
      }
      // loglikelihood
      double loglik_normalizing_const = - nk * (p+1.0)/2.0 * log(2.0*M_PI);
      loglikelihood += loglik_normalizing_const;
  } // end nk>0


  vec log_prior = logPriorDens( particles, varphi, prior );
  // Numerator of the unnormalized importance weights
  vec log_pi = log_prior + loglikelihood;
  return(log_pi);
}

Rcpp::List PMC::iter(uvec T,
                     int k,
                     umat N,
                     Rcpp::List particles,
                     mat log_dQ,
                     vec varphi,
                     Rcpp::List prior ) {

  uvec T_k = arma::find(T==k);
  mat Y_k = Y.rows(T_k);
  uvec C_k = C(T_k);

  int nk = Y_k.n_rows;

  // Proposal step (with permuted sweeps):

  List drawList;
  // cout << "start z" << endl;
  drawList = sampleZ(Y_k, C_k, particles);
  particles["z"] = Rcpp::as<mat>(drawList["z"]);
  log_dQ.col(1) = Rcpp::as<vec>(drawList["log_dq"]);
  // cout << "end z" << endl;

  // cout << "start S" << endl;
  drawList = sampleS( Y_k, C_k, N.col(k), particles, varphi);
  particles["S"] = Rcpp::as<vec>(drawList["S"]);
  log_dQ.col(2) = Rcpp::as<vec>(drawList["log_dq"]);
  // cout << "end S" << endl;

  uvec parPerm = randsamp(5,3,7);

  for(int ipar=0; ipar<5; ipar++) {
    switch ( parPerm(ipar) ) {
    case 3:
      // cout << "start xi" << endl;
      drawList = sampleXi( Y_k, C_k, N.col(k), particles);
      particles["xi"] = Rcpp::as<cube>(drawList["xi"]);
      log_dQ.col(3) = Rcpp::as<vec>(drawList["log_dq"]);
      // cout << "end xi" << endl;
      break;
    case 4:
      // cout << "start G" << endl;
      drawList = sampleG(Y_k, C_k, N.col(k), particles, prior);
      particles["G"] = Rcpp::as<mat>(drawList["G"]);
      log_dQ.col(4) = Rcpp::as<vec>(drawList["log_dq"]);
      // cout << "end G" << endl;
      break;
    case 5:
      // cout << "start psi" << endl;
      drawList = samplePsi( Y_k, C_k, N.col(k), particles, prior );
      particles["psi"] = Rcpp::as<mat>(drawList["psi"]);
      log_dQ.col(5) = Rcpp::as<vec>(drawList["log_dq"]);
      // cout << "end psi" << endl;
      break;
    case 6:
      // cout << "start xi0" << endl;
      drawList = sampleXi0( Y_k, N.col(k), particles);
      particles["xi0"] = Rcpp::as<mat>(drawList["xi0"]);
      log_dQ.col(6) = Rcpp::as<vec>(drawList["log_dq"]);
      // cout << "end xi0" << endl;
      break;
    case 7:
      // cout << "start E" << endl;
      drawList = sampleE( N.col(k), particles, prior);
      particles["E"] = Rcpp::as<mat>(drawList["E"]);
      log_dQ.col(7) = Rcpp::as<vec>(drawList["log_dq"]);
      // cout << "end E" << endl;
      break;
    }
  }

  vec iw;
  double log_py, perplexity;
  vec log_pitilde;
  if (nk>0) {
    log_pitilde = logPostDens(Y_k, C_k, N.col(k), particles, varphi, prior);
  } else {
    log_pitilde = logPriorDens( particles, varphi, prior );
  }
  vec log_q = sum(log_dQ, 1);
  // cout << "q " << log_q.t() << endl;
  vec log_iw = log_pitilde - log_q;
  // cout << "log iw raw " << log_iw.t() << endl;
  double cnst = max(log_iw);	// Constant needed for the computation of the marginal likelihood
  log_py = cnst + log(accu(exp(log_iw-cnst))) - log(num_particles);
  vec iw_bar = exp(log_iw);
  vec log_iw_b = log_iw - log(num_particles) - log_py;
  vec iw_b = exp(log_iw_b);  // unnormalized weights
  // cout << "log iw unnormalized " << iw_b.t() << endl;
  // cout << "iw unnormalized " << exp(iw_b.t()) << endl;
  iw = iw_b/accu(iw_b); // normalized weights
  // cout << "iw final " << iw.t() << endl;
  perplexity = exp(-accu(iw % log(iw))) / num_particles;


  // cout << "end copmute weights" << endl;


  // if (nk==0) cout << "start empty sample with replacement " << endl;

  uvec resamp(num_particles);
  for (int iN=0; iN<num_particles; iN++) {
    resamp(iN) = sampling(iw);
  }

  // if (nk==0) cout << "end empty sample with replacement " << endl;
  // cout << "iw " << iw.t() << endl;
  // NumericVector oo = Rcpp::NumericVector(iw.begin(), iw.end());
  // // cout << "oo " << oo << endl;
  // IntegerVector frame = seq(0, num_particles-1);
  // // cout << "frame " << frame << endl;
  // IntegerVector resampTemp = Rcpp::sample(frame, num_particles, true, oo);
  // // cout << "resamp " << resampTemp << endl;
  // uvec resamp = Rcpp::as<uvec>(resampTemp);
  // cout << "resamp " << resamp.t() << endl;
  uvec uresamp = unique(resamp);
  // cout << uJ << endl;
  int nResampled = uresamp.n_elem;
  // cout << "nresamp " << nResampled << endl;
  // cout << "end sample with replacement " << endl;

  // // Resampling step
  // cout << "start resample" << endl;
  cube MMM = Rcpp::as<cube>(particles["xi"]);
  cube MM = MMM;
  for (int i=0; i<num_particles; i++) {
    MM.slice(i) = MMM.slice(resamp(i));
  }
  particles["xi"] = MM;
  mat M;
  M = Rcpp::as<mat>(particles["G"]).rows(resamp);
  particles["G"] = M;
  // cout << M << endl;
  M = Rcpp::as<mat>(particles["psi"]).rows(resamp);
  particles["psi"] = M;
  // if (nk==0) cout << "start resamp xi0" << endl;
  M = Rcpp::as<mat>(particles["xi0"]).rows(resamp);
  particles["xi0"] = M;
  // if (nk==0) cout << "end resamp xi0" << endl;
  M = Rcpp::as<mat>(particles["E"]).rows(resamp);
  particles["E"] = M;

  // cout << "end resampling" << endl;

  return Rcpp::List::create(
    Rcpp::Named( "particles" ) = particles,
    Rcpp::Named( "log_py" ) = log_py,
    Rcpp::Named( "nResampled" ) = nResampled,
    Rcpp::Named( "perplexity" ) = perplexity  );
}

arma::uvec PMC::sampleT( arma::cube xi,
                         arma::cube Omega,
                         arma::mat alpha,
                         arma::mat logW )
{
  uvec T(n);
  int k;
  int j;
  uvec C_j;
  mat PT(n, K);
  PT.fill(0);
#pragma omp parallel for private(j, C_j)
  for (k=0; k<K; k++) {
    uvec k_idx(1);
    k_idx(0) = k;
    for (j=0; j<J; j++) {
      C_j = arma::find(C==j);
      PT.submat(C_j, k_idx) = dmsnArma( Y.rows(C_j),
                                        xi.slice(k).row(j),
                                        Omega.slice(k),
                                        alpha.col(k)) * exp(logW(j,k));
    }
  }

  NumericVector U = runif(n);
  double x;
  bool not_assigned;
  vec prob;
  vec probsum;
#pragma omp parallel for private(k, prob, probsum, x, not_assigned)
  for (int i=0; i<n; i++) {
    prob = PT.row(i).t();
    probsum = cumsum(prob);
    x = U(i) * sum(prob);
    not_assigned = true;
    for (k=0; (k<K) && not_assigned; k++) {
      if (x <= probsum(k)) {
        T(i) = k;
        not_assigned = false;
      }
    }
  }

  return(T);
}

//--------------------------------------------------------

Rcpp::List PMC::get_chain()
{
  return Rcpp::List::create(
    Rcpp::Named( "t" ) = saveT,
    Rcpp::Named( "z" ) = saveZ,
    Rcpp::Named( "W" ) = saveW,
    Rcpp::Named( "xi" ) = saveXi,
    Rcpp::Named( "xi0" ) = saveXi0,
    Rcpp::Named( "psi" ) = savePsi,
    Rcpp::Named( "G" ) = saveG,
    Rcpp::Named( "E" ) = saveE,
    Rcpp::Named( "S" ) = saveS,
    Rcpp::Named( "varphi" ) = saveVarphi,
    Rcpp::Named( "a0" ) = saveA0,
    Rcpp::Named( "log_py" ) = saveLog_py,
    Rcpp::Named( "perplexity" ) = savePerplexity,
    Rcpp::Named( "nResampled" ) = saveNResampled,
    Rcpp::Named( "Omega" ) = saveOmega,
    Rcpp::Named( "alpha" ) = saveAlpha
  );
};
