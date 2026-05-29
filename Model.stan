functions {
  matrix GP_cor_cholesky(matrix dist_mat_sq, real rho) {
    int N = rows(dist_mat_sq);
    matrix[N, N] K = exp(-dist_mat_sq / (2.0 * square(rho)));
    return cholesky_decompose(add_diag(K, 1e-9));
  }

  // UPDATED: Used array[] int syntax for function arguments
  vector simulate_abundance(int S, int T_total, int T_obs,
                            vector thermal_niche, array[] int temp_idx, 
                            vector phi, vector gamma_smooth, vector t_scaled, 
                            vector eps_year, real r, real N0_proportion, 
                            vector w_norm, array[] int from_idx, array[] int to_idx, int E_disp) {
    
    int N_state_total = S * T_obs;
    vector[N_state_total] N_state_vec;
    real R0 = exp(r);
    real R0m1 = R0 - 1.0;
    
    vector[S] N_curr;
    
    // Initial state (Year 1)
    {
      vector[S] logN = thermal_niche[temp_idx[1:S]] + phi + (gamma_smooth * t_scaled[1]) + eps_year[1];
      N_curr = exp(logN) * N0_proportion;
    }
    
    // IDE Loop
    for (t in 2:T_total) {
      vector[S] logK = thermal_niche[temp_idx[((t - 1) * S + 1) : (t * S)]] + phi + (gamma_smooth * t_scaled[t]);
      vector[S] K = exp(logK);
      
      // Beverton-Holt growth
      N_curr = (R0 * N_curr) ./ (1.0 + (R0m1 ./ (K + 1e-9)) .* N_curr);
      
      // Dispersal
      vector[S] N_new = rep_vector(0.0, S);
      for (e in 1:E_disp) N_new[to_idx[e]] += w_norm[e] * N_curr[from_idx[e]];
      
      N_curr = N_new * exp(eps_year[t]);
      
      // Store state if within observation period
      if (t > T_total - T_obs) {
        int obs_t = t - (T_total - T_obs);
        N_state_vec[((obs_t - 1) * S + 1):(obs_t * S)] = N_curr;
      }
    }
    return N_state_vec;
  }
}

data {
  int<lower=1> S;
  int<lower=1> T_total;
  int<lower=1> T_obs;
  int<lower=1> N_obs_total;
  int<lower=1> N_state_total;
  vector[T_total] t_scaled;
  array[N_obs_total] int<lower=0> Y_vec;
  array[N_obs_total] int<lower=0> K_visits_vec;
  array[N_obs_total] int<lower=1, upper=N_state_total> map_state_idx;
  int<lower=1> E_disp;
  vector[E_disp] dists;
  array[E_disp] int from_idx;
  array[E_disp] int to_idx;
  int<lower=1> N_spatial_bf;
  matrix[S, N_spatial_bf] spatial_bf;
  matrix[N_spatial_bf, N_spatial_bf] dist_mat_anchors_sq;
  int<lower=1> N_thermal_gradient;
  vector[N_thermal_gradient] thermal_gradient;
  array[S * T_total] int<lower=1, upper=N_thermal_gradient> temp_idx;
}

transformed data {
  int N_pos = 0;
  int N_zero = 0;
  for (i in 1:N_obs_total) {
    if (Y_vec[i] > 0) N_pos += 1;
    else N_zero += 1;
  }
  
  array[N_pos] int pos_state_idx;
  array[N_pos] int pos_Kvis;
  array[N_pos] int pos_Y;
  array[N_zero] int zero_state_idx;
  array[N_zero] int zero_Kvis;

  int p_ptr = 1;
  int z_ptr = 1;
  for (i in 1:N_obs_total) {
    if (Y_vec[i] > 0) {
      pos_state_idx[p_ptr] = map_state_idx[i];
      pos_Kvis[p_ptr] = K_visits_vec[i];
      pos_Y[p_ptr] = Y_vec[i];
      p_ptr += 1;
    } else {
      zero_state_idx[z_ptr] = map_state_idx[i];
      zero_Kvis[z_ptr] = K_visits_vec[i];
      z_ptr += 1;
    }
  }
}

parameters {
  real<lower=0> r;
  real log_alpha;
  real<lower=0, upper=1> N0_proportion;
  real<lower=0> sigma_year;
  vector[T_total] eps_year_raw;
  real<lower=0> rho_phi;
  real<lower=0> rho_gamma;
  real<lower=0> sigma_phi;
  vector[N_spatial_bf] phi_eta;
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
  eps_year -= mean(eps_year);

  matrix[N_spatial_bf, N_spatial_bf] L_phi = GP_cor_cholesky(dist_mat_anchors_sq, rho_phi);
  matrix[N_spatial_bf, N_spatial_bf] L_gamma = GP_cor_cholesky(dist_mat_anchors_sq, rho_gamma);

  vector[S] phi = sigma_phi * (spatial_bf * (L_phi * phi_eta));
  phi -= mean(phi);

  vector[S] gamma_smooth = spatial_bf * (L_gamma * gamma_eta);
  gamma_smooth -= mean(gamma_smooth); 

  vector[N_thermal_gradient] thermal_niche;
  for (k in 1:N_thermal_gradient) {
    real x = thermal_gradient[k];
    thermal_niche[k] = thermal_max - log1p_exp(-slope_L * (x - T_nodes[1])) - log1p_exp(-slope_R * (T_nodes[2] - x));
  }
}

model {
  r ~ normal(0.3, 0.3);
  log_alpha ~ normal(2.5, 1);
  N0_proportion ~ beta(2, 2);
  sigma_year ~ student_t(3, 0, 0.2);
  eps_year_raw ~ std_normal();
  rho_phi ~ lognormal(log(100), 0.5);
  rho_gamma ~ lognormal(log(150), 0.5);
  sigma_phi ~ normal(0, 0.1);
  phi_eta ~ std_normal();
  sigma_gamma ~ normal(0, 0.1);
  gamma_eta ~ std_normal();
  logit_p ~ normal(-2, 1);
  thermal_max ~ normal(0, 3);
  T_nodes ~ normal(10, 8);
  slope_L ~ lognormal(0, 0.5);
  slope_R ~ lognormal(0, 0.5);

  vector[E_disp] w = exp(-dists / alpha);
  vector[S] out_sum = rep_vector(1e-12, S);
  for (e in 1:E_disp) out_sum[from_idx[e]] += w[e];
  vector[E_disp] w_norm = w ./ out_sum[from_idx];

  vector[N_state_total] N_state_vec = simulate_abundance(S, T_total, T_obs, 
                                        thermal_niche, temp_idx, phi, gamma_smooth, 
                                        t_scaled, eps_year, r, N0_proportion, 
                                        w_norm, from_idx, to_idx, E_disp);

  vector[N_state_total] log_psi = log1m_exp(-(N_state_vec + 1e-12));
  
  target += sum(log_psi[pos_state_idx]) + binomial_logit_lpmf(pos_Y | pos_Kvis, logit_p);

  for (i in 1:N_zero) {
    int s = zero_state_idx[i];
    target += log_sum_exp(-N_state_vec[s], log_psi[s] + binomial_logit_lpmf(0 | zero_Kvis[i], logit_p));
  }
}

generated quantities {
  vector[E_disp] w = exp(-dists / alpha);
  vector[S] out_sum = rep_vector(1e-12, S);
  for (e in 1:E_disp) out_sum[from_idx[e]] += w[e];
  vector[E_disp] w_norm = w ./ out_sum[from_idx];

  vector[N_state_total] N_state = simulate_abundance(S, T_total, T_obs, 
                                        thermal_niche, temp_idx, phi, gamma_smooth, 
                                        t_scaled, eps_year, r, N0_proportion, 
                                        w_norm, from_idx, to_idx, E_disp);
}
