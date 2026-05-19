functions {
  matrix GP_cor_cholesky(matrix dist_mat, real rho) {
    int N = rows(dist_mat);
    matrix[N, N] K = exp(-square(dist_mat) / (2 * square(rho)));
    return cholesky_decompose(add_diag(K, 1e-9));
  }

  vector log_sum_exp_vector(vector a, vector b) {
    vector[rows(a)] res;
    for (i in 1:rows(a)) {
      res[i] = log_sum_exp(a[i], b[i]);
    }
    return res;
  }
}

data {
  int<lower=1> S;
  int<lower=1> T_total;
  int<lower=1> T_obs;
  int<lower=1> N_obs_total;
  array[N_obs_total] int<lower=0> Y_vec;
  array[N_obs_total] int<lower=0> K_visits_vec;
  int<lower=1> N_state_total; 
  array[N_obs_total] int<lower=1, upper=N_state_total> map_state_idx;

  int<lower=0> N_seen;
  int<lower=0> N_zero;
  array[N_seen] int idx_seen;
  array[N_zero] int idx_zero;

  int<lower=1> E_disp;
  vector[E_disp] dists;
  array[E_disp] int to_idx;
  array[S+1] int row_ptr;
  array[E_disp] int row_ids; // NEW
  
  int<lower=1> N_spatial_bf;
  matrix[S, N_spatial_bf] spatial_bf;
  matrix[N_spatial_bf, N_spatial_bf] dist_mat_anchors; 
  
  int<lower=1> N_thermal_gradient;
  vector[N_thermal_gradient] thermal_gradient;
  array[S * T_total] int<lower=1, upper=N_thermal_gradient> temp_idx;
}

transformed data {
  vector[T_total] time_scaled;
  for (t in 1:T_total) {
    time_scaled[t] = (t - T_total / 2.0) / T_total;
  }
}

parameters {
  real<lower=0> r;
  real log_alpha;
  real<lower=0, upper=1> N0_proportion;
  real<lower=0> sigma_year;
  vector[T_total] eps_year_raw;
  
  real<lower=0> rho_phi;
  real<lower=0> sigma_phi;
  vector[N_spatial_bf] phi_eta;
  
  real<lower=0> rho_gamma; // Ensure you have rho_gamma defined!
  real<lower=0> sigma_gamma;        
  vector[N_spatial_bf] gamma_eta;   

  real logit_p;    
  
  real thermal_max;             
  ordered[2] T_nodes;           
  real<lower=0> slope_L;        
  real<lower=0> slope_R;        
}

transformed parameters {
  real alpha = exp(log_alpha);
  vector[T_total] eps_year = sigma_year * eps_year_raw;
  
  // 1. Centered Habitat Field
  vector[S] phi;
  {
    matrix[N_spatial_bf, N_spatial_bf] L_spatial = GP_cor_cholesky(dist_mat_anchors, rho_phi);
    vector[S] phi_raw = spatial_bf * (sigma_phi * (L_spatial * phi_eta));
    phi = phi_raw - mean(phi_raw);
  }

  // 2. Centered Redistribution Trend
  vector[S] gamma;
  {
    matrix[N_spatial_bf, N_spatial_bf] L_gamma = GP_cor_cholesky(dist_mat_anchors, rho_gamma);
    vector[S] gamma_raw = spatial_bf * (sigma_gamma * (L_gamma * gamma_eta));
    gamma = gamma_raw - mean(gamma_raw);
  }
  
  vector[N_state_total] N_state_vec; 
  {
    matrix[S, 2] N_buffer;
    real R0 = exp(r); 

    // 3. Uncentered Thermal Niche
    vector[N_thermal_gradient] static_niche = thermal_max 
           + log_inv_logit(slope_L * (thermal_gradient - T_nodes[1])) 
           + log_inv_logit(slope_R * (T_nodes[2] - thermal_gradient));

    // Precompute Niche Lookup
    matrix[T_total, S] niche_mat;
    for (t in 1:T_total) {
      int start = (t-1)*S + 1;
      niche_mat[t] = to_row_vector(static_niche[temp_idx[start:(start+S-1)]]);
    }

    vector[E_disp] w_norm;
    {
       vector[E_disp] w = exp(-dists / alpha);
       vector[S] out_sums = rep_vector(1e-9, S);
       
       // [FIX 1]: Row-stochastic normalization!
       for (e in 1:E_disp) {
         out_sums[row_ids[e]] += w[e];
       }
       w_norm = w ./ out_sums[row_ids];
    }

    // Initialize Year 1 (Reverting to simple bounds)
    vector[S] log_init = to_vector(niche_mat[1]') + phi + gamma * time_scaled[1] + eps_year[1];
    N_buffer[,1] = fmin(fmax(exp(log_init), 1e-6), 1e5) * N0_proportion;

    for (t in 2:T_total) {
      vector[S] log_K_raw = to_vector(niche_mat[t]') + phi + gamma * time_scaled[t];
      vector[S] K = fmin(fmax(exp(log_K_raw), 1e-6), 1e5);
      
      vector[S] growth = (R0 * N_buffer[,1]) ./ (1.0 + ((R0 - 1.0) ./ (K + 1e-9)) .* N_buffer[,1]);
      
      // CSR Dispersal
      vector[S] dispersed_raw = csr_matrix_times_vector(S, S, w_norm, to_idx, row_ptr, growth);
      
      // Hard cap handles positivity directly (Massive Speedup)
      N_buffer[,2] = fmin(fmax(dispersed_raw .* exp(eps_year[t]), 1e-6), 1e5);

      if (t > T_total - T_obs) {
        int obs_t = t - (T_total - T_obs);
        N_state_vec[((obs_t-1)*S + 1):(obs_t*S)] = N_buffer[,2]; 
      }
      N_buffer[,1] = N_buffer[,2]; 
    }
  }
}

model {
  // GP Lengthscales
  rho_phi ~ lognormal(log(200), 0.5);   // Expect habitat smoothness over ~200km
  rho_gamma ~ lognormal(log(500), 0.4); // Expect temporal trends over ~500km
  
  // GP Magnitudes (Strictly one prior per parameter)
  sigma_phi ~ student_t(3, 0, 0.25); 
  sigma_gamma ~ normal(0, 0.1); 
  
  // Standard Normal for Latent Fields
  phi_eta ~ std_normal();
  gamma_eta ~ std_normal();
  
  // Safer Niche Priors
  thermal_max ~ normal(0, 3); 
  T_nodes[1] ~ normal(5, 5);  
  T_nodes[2] ~ normal(15, 5); 
  slope_L ~ lognormal(0, 0.5); 
  slope_R ~ lognormal(0, 0.5); 
  
  // Ecological & Temporal Priors
  r ~ normal(0, 0.3); 
  log_alpha ~ normal(2.7, 1); 
  sigma_year ~ student_t(3, 0, 0.25);
  eps_year_raw ~ std_normal(); 
  N0_proportion ~ beta(2, 2); 
  
  // Realistic Detection Prior
  logit_p ~ normal(-2.0, 1.0);     

  // Likelihood
  vector[N_state_total] log_psi = log1m_exp(-N_state_vec); 
  vector[N_state_total] log1m_psi = -N_state_vec; 

  target += log_psi[map_state_idx[idx_seen]];
  target += binomial_logit_lpmf(Y_vec[idx_seen] | K_visits_vec[idx_seen], logit_p);

  real log_prob_miss = log1m_inv_logit(logit_p);
  target += sum(log_sum_exp_vector(
      log1m_psi[idx_zero], 
      log_psi[idx_zero] + to_vector(K_visits_vec[idx_zero]) * log_prob_miss
  ));
}

