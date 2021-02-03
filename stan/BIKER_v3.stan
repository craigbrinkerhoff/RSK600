//Craig BEE stan model for remotely sensing gas transfer velocity

functions {
  // Conversion from array to vector takes place by row.
  // Nested elements are appended in sequence.

  // Convert an array to a vector based on a binary matrix
  // indicating non-missing data
  vector ragged_vec(vector[] x, int[,] bin) {
    vector[num_elements(x)] out;
    int ind;

    ind = 1;
    for (i in 1:size(x)) {
      for (t in 1:num_elements(x[1])) {
        if (bin[i, t] == 1) {
          out[ind] = x[i, t];
          ind += 1;
        }
      }
    }
    // print(out);
    return(out[1:(ind - 1)]);
  }

  // Repeat elements of a "row" vector to match with 2-D array vectorization
  vector ragged_row(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;

    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[t];
        }
      }
    }
    return(out[1:ind]);
  }

  // Repeat elements of a "column" vector to match with 2-D array vectorization
  vector ragged_col(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;

    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[i];
        }
      }
    }
    return(out[1:ind]);
  }

  // indices of vectorized bin1 that are also in vectorized bin2
  int[] commoninds(int[,] bin1, int[,] bin2) {
    int out[num_elements(bin1)];
    int vecinds[size(bin1), num_elements(bin1[1])];
    int ctr;
    int ind;

    ctr = 0;
    for (i in 1:size(bin1)) {
      for (t in 1:num_elements(bin1[1])) {
        if (bin1[i, t] == 1) {
          ctr += 1;
          vecinds[i, t] = ctr;
        }
      }
    }

    ind = 0;
    for (i in 1:size(vecinds)) {
      for (t in 1:size(vecinds[1])) {
        if (bin2[i, t] == 1) {
          ind += 1;
          out[ind] = vecinds[i, t];
        }
      }
    }
    return(out[1:ind]);
  }
}

data {

  // Options
  int<lower=0, upper=1> inc_m; // Include Manning? 0=no; 1=yes;

  // Dimensions
  int<lower=0> nx; // number of reaches
  int<lower=0> nt; // number of times
  int<lower=0> ntot_man; // total number of non-missing Manning observations

  // // Missing data
  int<lower=0,upper=1> hasdat_man[nx, nt]; // matrix of 0 (missing), 1 (not missing)

  // *Actual* data
  vector[nt] Wobs[nx]; // measured widths, including placeholders for missing
  vector[nt] Sobs[nx]; // measured slopes
  vector[nt] dAobs[nx]; // measured partial area
  vector[nx] dA_shift; // adjustment from min to median

  // Hard bounds on parameters
  real lowerbound_A0; // These must be scalars, unfortunately.
  real upperbound_A0;
  real lowerbound_logn;
  real upperbound_logn;
  real upperbound_logk600;
  real lowerbound_logk600;

  // *Known* likelihood parameters
  vector<lower=0>[nt] sigma_post[nx]; // Manning error standard deviation

  // Hyperparameters
  vector[nt] logk600_hat; // prior mean on logF g/m2*dy
  real logA0_hat[nx]; //space-varying A0 prior m2
  real logn_hat[nx]; //space varying n prior m2

  vector<lower=0>[nt] logk600_sd;
  real<lower=0> logA0_sd[nx];
  real<lower=0> logn_sd[nx];
}



transformed data {
  // Transformed data are *vectors*, not arrays. This to allow ragged structure

  vector[nt] dApos_array[nx];

  vector[ntot_man] Wobsvec_man;
  vector[ntot_man] Sobsvec_man;

  vector[ntot_man] logWobs_man;
  vector[ntot_man] logSobs_man;
  vector[ntot_man] dApos_obs;
  vector[ntot_man] sigmavec_man;

  int ntot_w; // how many widths in likelihood: equal to ntot_man unless inc_a
  ntot_w = ntot_man;

  for (i in 1:nx) {
    dApos_array[i] = dAobs[i] - min(dAobs[i]); // make all dA positive
  }

  // convert pseudo-ragged arrays to vectors
  Wobsvec_man = ragged_vec(Wobs, hasdat_man);
  Sobsvec_man = ragged_vec(Sobs, hasdat_man);
  dApos_obs = ragged_vec(dApos_array, hasdat_man);

  logWobs_man = log(Wobsvec_man);
  logSobs_man = log(Sobsvec_man);

  sigmavec_man = ragged_vec(sigma_post, hasdat_man);
}

parameters {
  vector<lower=lowerbound_logk600,upper=upperbound_logk600>[nt] logk600;
  vector<lower=lowerbound_logn,upper=upperbound_logn>[nx] logn[inc_m]; //for reach-defined n
  vector<lower=lowerbound_A0,upper=upperbound_A0>[nx] A0[inc_m];
}


transformed parameters {

  vector[ntot_man] man_lhs[inc_m]; // LHS for Manning likelihood
  vector[ntot_man] logA_man[inc_m]; // log area for Manning's equation
  vector[ntot_man] logN_man[inc_m]; //log Manning's n for Manning's equation
  vector[ntot_man] logk600_man[inc_m]; // location-repeated logk600
  vector[ntot_man] man_rhs[inc_m]; // RHS for Manning likelihood

  // Manning params
  if (inc_m) {
    logA_man[1] = log(ragged_col(A0[1], hasdat_man) + dApos_obs);
    logN_man[1] = ragged_col(logn[1], hasdat_man);
    logk600_man[1] = ragged_row(logk600, hasdat_man);

    //Rule-based (using average obs width) model for k600. Basically, different regression parameters depending upon river size
    if(mean(logWobs_man) < log(10) && mean(Sobsvec_man) < 0.05){
      man_lhs[1] = 0.4325421*logWobs_man - 0.9732197*logSobs_man - log(111.58121) - 0.6488131*log(9.8);
      man_rhs[1] = 0.4325421*(logA_man[1]) - 0.6488131*logN_man[1] - logk600_man[1];
    }
    if(mean(logWobs_man) < log(10) && mean(Sobsvec_man) >= 0.05){
      man_lhs[1] = 0.8773377*logWobs_man - 1.97401*logSobs_man - log(792.63149) - 1.3160065*log(9.8);
      man_rhs[1] = 0.8773377*(logA_man[1]) - 1.3160065*logN_man[1] - logk600_man[1];
    }
    if(mean(logWobs_man) < log(50) && mean(logWobs_man) >= log(10)){
      man_lhs[1] = 0.4416236*logWobs_man - 0.9936531*logSobs_man - log(109.04977) - 0.6624354*log(9.8);
      man_rhs[1] = 0.4416236*(logA_man[1]) - 0.6624354*logN_man[1] - logk600_man[1];
    }
    if(mean(logWobs_man) < log(100) && mean(logWobs_man) >= log(50)){
      man_lhs[1] = 0.2910743*logWobs_man - 0.6549171*logSobs_man - log(31.84344) - 0.4366114*log(9.8);
      man_rhs[1] = 0.2910743*(logA_man[1]) - 0.4366114*logN_man[1] - logk600_man[1];
    }
    if(mean(logWobs_man) >= log(100)){
      man_lhs[1] = 0.1816556*logWobs_man - 0.4087251*logSobs_man - log(14.16939) - 0.2724834*log(9.8);
      man_rhs[1] = 0.1816556*(logA_man[1]) - 0.2724834*logN_man[1] - logk600_man[1];
    }
  }
}

model {
  // Priors
  if (inc_m) {
    A0[1] + dA_shift[1] ~ lognormal(logA0_hat, logA0_sd);
    logn[1] ~ normal(logn_hat, logn_sd);
    logk600[1] ~ normal(logk600_hat, logk600_sd);
}
  //Run actual model
  if (inc_m){
      man_lhs[1] ~ normal(man_rhs[1], sigmavec_man);
  }
}
