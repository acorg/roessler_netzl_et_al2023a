
functions {
  
  real normal_int_censored_likelihood(real lower_lim, real upper_lim, real mu, real sigma) {
    
    real result;
    
    if (is_inf(lower_lim) && is_inf(upper_lim)) {
      
      result = 0;
      
    } else if (lower_lim == upper_lim) {
      
      result = normal_lpdf(
        lower_lim | mu, sigma
      );
      
    } else if (!is_inf(lower_lim) && !is_inf(upper_lim)) {
  
      result = log_diff_exp(
        normal_lcdf(upper_lim | mu, sigma),
        normal_lcdf(lower_lim | mu, sigma)
      );
      
    } else if (!is_inf(lower_lim) && is_inf(upper_lim)) {
      
      result = normal_lccdf(
        lower_lim | mu, sigma
      );
      
    } else if (is_inf(lower_lim) && !is_inf(upper_lim)) {
      
      result = normal_lcdf(
        upper_lim | mu, sigma
      );
      
    }
    
    return result;
    
  }
  
  real normal_lub_rng(real mu, real sigma, real lower_lim, real upper_lim) {
    
    if (lower_lim == upper_lim) {
      
      return lower_lim;
      
    } else {
      
      real p_lower_lim = normal_cdf(lower_lim, mu, sigma);
      real p_upper_lim = normal_cdf(upper_lim, mu, sigma);
      real u = uniform_rng(p_lower_lim, p_upper_lim);
      real y = mu + sigma * inv_Phi(u);
      return y;
      
    }
  }
  
}

data {
  int<lower=1> N;
  int<lower=1> N_ags;
  int<lower=1> N_srs;
  int<lower=1> N_sr_groups;
  int<lower=1> N_datasets;
  vector[N] upper_lims;
  vector[N] lower_lims;
  array[N] int ags;
  array[N] int srs;
  array[N] int sr_groups;
  array[N] int datasets;
  real ag_mean_prior_mean;
  real<lower=0> ag_mean_prior_sigma;
  real<lower=0> sigma_prior_alpha;
  real<lower=0> sigma_prior_beta;
  real<lower=0> sigma_prior_sr_effect;
  real<lower=0> sigma_prior_dataset_magnitude;
  vector[N_datasets] dataset_magnitude_mean;
}

parameters {
  matrix[N_ags, N_sr_groups] ag_means;
  vector[N_srs] sr_effects;
  vector[N_datasets] dataset_magnitude;
  matrix<lower=0>[N_datasets, N_sr_groups] sigma_error;
}

model {
  
  // Calculate likelihood of priors
  sr_effects ~ normal(0, sigma_prior_sr_effect);
  
  for (i in 1:N_datasets) {
    dataset_magnitude[i] ~ normal(dataset_magnitude_mean[i], sigma_prior_dataset_magnitude);
  }
  
  
  for (i in 1:N_sr_groups) {
    sigma_error[,i] ~ inv_gamma(sigma_prior_alpha, sigma_prior_beta);
    ag_means[,i] ~ normal(ag_mean_prior_mean, ag_mean_prior_sigma);
  }
  
  
  // Work out likelihood of each titer
  for (i in 1:N) {
    
    // Set variables
    int ag = ags[i];
    int sr = srs[i];
    int sr_group = sr_groups[i];
    int dataset = datasets[i];
    
    // Add in mixed effect for each antigen and serum
    real logtiter = ag_means[ag, sr_group] + 
      dataset_magnitude[dataset] + 
      sr_effects[sr]; 
      
    
    // Work out titer likelihood
    target += normal_int_censored_likelihood(
      lower_lims[i],
      upper_lims[i],
      logtiter, 
      sigma_error[dataset, sr_group]
    );
    
  }
  
}

generated quantities {
  
  vector[N] imputed_logtiters;
  
  // Work out likelihood of each titer
  for (i in 1:N) {
    
    // Set variables
    int ag = ags[i];
    int sr = srs[i];
    int sr_group = sr_groups[i];
    int dataset = datasets[i];
    
    // Add in mixed effect for each antigen and serum
     real logtiter = ag_means[ag, sr_group] + 
      dataset_magnitude[dataset] + 
      sr_effects[sr]; 

      // Add in mixed effect for each antigen and serum
    imputed_logtiters[i] = normal_lub_rng(
      logtiter,
      sigma_error[dataset, sr_group],
      lower_lims[i],
      upper_lims[i]
      );
  
   
    
    
  }
  
}