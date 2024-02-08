#*******************************************************************************
#*******************************************************************************
#* The functions used for the paper "Alone, together: on the benefits of 
#* Bayesian borrowing in a meta-analytic setting", to appear on Pharmaceutical
#* Statistics (2023).
#* 
#* All code was written by Ofir Harari
#*******************************************************************************
#*******************************************************************************




#*******************************************************************************
#* Generating random (aggregate-level) data
#* 
#* Inputs: 
#*   K - the number of studies
#*   p_ctrl - control event rate
#*   p_trt - treatment event rate
#*   sigma2_ctrl - variance for the control group study-to-study perturbation
#*   sigma2_trt - variance for the treatment group study-to-study perturbation
#*   n_ctrl - control group number of patients (per study)
#*   n_trt - treatment group number of patients (per study)
#*   
#* Output: a list consisting of true and estimated treatment effects and their 
#*         squared standard errors.
#*******************************************************************************
random_data = function(K, p_ctrl, p_trt, sigma2_ctrl, 
                       sigma2_trt, n_ctrl, n_trt){
  p_ctrl_vec = LaplacesDemon::invlogit(rnorm(K, logit(p_ctrl), sigma2_ctrl))
  p_trt_vec = LaplacesDemon::invlogit(rnorm(K, logit(p_trt), sigma2_trt))
  theta_vec  = log(p_trt_vec/p_ctrl_vec)
  
  n_ctrl_vec = rpois(K, 75)
  n_trt_vec = rpois(K, 75)
  y_ctrl_vec = rbinom(K, n_ctrl_vec, p_ctrl_vec)
  y_trt_vec = rbinom(K, n_trt_vec, p_trt_vec)
  theta_hat_vec  = log(y_trt_vec*n_ctrl_vec/y_ctrl_vec/n_trt_vec)
  
  sigma2_vec = 1/y_ctrl_vec + 1/y_trt_vec - 1/n_ctrl_vec - 1/n_trt_vec
  
  return(list(theta_vec = theta_vec, 
              theta_hat_vec = theta_hat_vec,
              sigma2_vec = sigma2_vec))
}



#*******************************************************************************
#* Displaying the random data that was generated
#* 
#* Inputs: 
#*   data - the output of the random_data() function
#*   
#* Output: a table
#*******************************************************************************
display_simulated_data = function(data){
  theta_vec = data$theta_vec
  theta_hat_vec = data$theta_hat_vec
  sigma2_vec = data$sigma2_vec
  
  cbind(theta_vec, theta_hat_vec, sigma2_vec) %>% 
    as.data.frame %>% 
    dplyr::rename('$\\log \\mathrm{RR}$' = 'theta_vec',
                  '$\\log \\widehat{\\mathrm{RR}}$' = 'theta_hat_vec',
                  '$\\mathrm{SE}^2$' = 'sigma2_vec') %>% 
    apply(., 2, round, digits = 3) %>% 
    kable(align = 'c') %>%
    kable_styling(full_width = F, font_size = 14,
                  bootstrap_options = c("striped")) %>%
    row_spec(0, font_size = 14) 
}



#*******************************************************************************
#* Calculating Bayesian borrowing fractions
#* Inputs: 
#*   sigma2_vec - the squared standard errors of the individual studies
#*   nu2 - the study-to-study variance
#*   i - the index of the target study to be updated by borrowing information 
#*       from other studies
#*   
#* Output: a vector of Bayesian borrowing fractions
#*******************************************************************************
weights = function(sigma2_vec, nu2, i){
  v = 1/(sigma2_vec + nu2)
  u = sigma2_vec[i]*v
  w = sigma2_vec/sigma2_vec[i]*u/(u[i] + nu2*sum(v))
  w[i] = 1
  
  w
}



#*******************************************************************************
#* Plotting the borrowing fraction per study vs. the study-to-study variance
#* 
#* Inputs: 
#*   sigma2_vec - the squared standard errors of the individual studies
#*   nu2_min - minimum nu2 value for evaluation
#*   nu2_max - maximum nu2 value for evaluation
#*   nu2_delta - distance between grid points
#*   trial_number - the index of the target study to be updated by borrowing 
#*                  information from other studies
#*   target_nu2 - the value of nu2 to be highlighted in the plot
#*   
#* Output: a ggplot
#*******************************************************************************
weight_vs_nu2_plot = function(sigma2_vec, nu2_min, nu2_max, 
                              nu2_delta, trial_number, target_nu2){
  nu2_vec = seq(nu2_min, nu2_max, by = nu2_delta)
  w = t(sapply(nu2_vec, weights, sigma2_vec = sigma2_vec, i = trial_number))
  colnames(w) = paste0(1:ncol(w))
  w = melt(w)[,-1]
  names(w) = c('Trial', 'Weight')
  w = w %>% 
    filter(Trial != trial_number) %>% 
    mutate(Trial = paste0("Trial ", Trial))
  w$nu2 = rep(nu2_vec, length(sigma2_vec) - 1)
  
  df_line = w %>% 
    filter(round(nu2, 4) == target_nu2) %>% 
    mutate(Text1 = paste0('nu^2 == ', nu2),
           Text2 = paste0('omega == ', sprintf('%.2f', Weight)))
  
  p = ggplot(w, aes(nu2, Weight)) + 
    geom_line(col = 'royal blue', linewidth = 1) + 
    facet_rep_grid(.~ Trial, repeat.tick.labels = 'all') +
    theme(panel.background = element_blank(),
          legend.position = 'top',
          panel.grid.major = element_line(colour = "grey92"),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=20, face="bold"),
          axis.title.y = element_text(size=18, face="bold"),
          axis.text.y = element_text(size=14),
          strip.text.x = element_text(size = 12, face="bold"),
          axis.line=element_line()) + 
    scale_x_continuous(expand = expansion(c(0,0.15))) + 
    geom_segment(df_line, 
                 mapping = aes(x = min(w$nu2), xend = nu2,
                               y = Weight, yend = Weight),
                 linetype = 'dashed') + 
    geom_segment(df_line, 
                 mapping = aes(x = nu2, xend = nu2,
                               y = Weight, yend = -Inf),
                 linetype = 'dashed') + 
    geom_point(df_line, mapping = aes(x = nu2, y = Weight),
               shape = 21, size = 3, fill = 'white', 
               stroke = 2, col = 'dark red') + 
    geom_text(df_line, mapping = aes(min(w$nu2) + .05, max(w$Weight)*.95, 
                                     label = Text1),
              size = 6, parse = T) + 
    geom_text(df_line, mapping = aes(min(w$nu2) + .05, max(w$Weight)*.85, 
                                     label = Text2),
              size = 6, parse = T) + 
    xlab(expression(nu^2)) + 
    ylab(expression('Borrowing fraction '~(omega)))
  
  return(p)
}




#*******************************************************************************
#* Posterior mean coefficients (of the individual treatment effects)
#* 
#* Inputs: 
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   nu2 - the study-to-study variance
#*   i - the index of the target study to be updated by borrowing information 
#*       from other studies
#*   
#* Output: a vector of coefficients (summing to 1)
#*******************************************************************************
a_coeffs = function(theta_hat_vec, 
                    sigma2_vec, nu2, i){
  k = length(sigma2_vec)
  v = 1/(sigma2_vec + nu2)
  sum_v = sum(v)
  a = sigma2_vec[i]*v[i]*v
  a[i] = a[i] + nu2*v[i]*sum_v
  
  return(a/sum_v)
}




#*******************************************************************************
#* The updated study-specific treatment effect estimate after Bayesian borrowing
#* 
#* Inputs: 
#*   theta_hat_vec - a vector of treatment effect estimated before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   nu2 - the study-to-study variance
#*   i - the index of the target study to be updated by borrowing information 
#*       from other studies
#*   
#* Output: an updated treatment effect estimate
#*******************************************************************************
adjusted_estimate = function(theta_hat_vec, 
                             sigma2_vec, nu2, i){
  k = length(sigma2_vec)
  v = 1/(sigma2_vec + nu2)
  sum_v = sum(v)
  a = sigma2_vec[i]*v[i]*v
  a[i] = a[i] + nu2*v[i]*sum_v
  
  return(a%*%theta_hat_vec/sum_v)
}



#*******************************************************************************
#* The updated study-specific treatment effect posterior variance after Bayesian 
#* borrowing
#* 
#* Inputs: 
#*   sigma2_vec - the squared standard errors of the individual studies
#*   nu2 - the study-to-study variance
#*   i - the index of the target study to be updated by borrowing information 
#*       from other studies
#*   
#* Output: posterior treatment effect variance
#*******************************************************************************
adjusted_variance = function(sigma2_vec, nu2, i){
  v = 1/(sigma2_vec + nu2)
  sum_v = sum(v)
  term = sigma2_vec[i]*v[i]
  
  term*(term + nu2*sum_v)/sum_v
}



#*******************************************************************************
#* Study-specific Bayesian update of all treatment effects jointly, for a fixed
#* study-to-study variance
#* 
#* Inputs: 
#*   theta_hat_vec - a vector of treatment effect estimated before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   alpha - 1 - the coverage level of the Bayesian credible intervals
#*   nu2 - the study-to-study variance
#*   
#* Output: data frame containing treatment effect estimates and uncertainty 
#*         intervals, before and after Bayesian borrowing
#*******************************************************************************
output_fixed_nu = function(theta_hat_vec, sigma2_vec, alpha, nu2){
  k = length(theta_hat_vec)
  
  theta_hat_vec_adj = sapply(1:k, adjusted_estimate,
                             theta_hat_vec = theta_hat_vec, 
                             sigma2_vec = sigma2_vec, nu2 = nu2)
  se = sqrt(sigma2_vec)
  
  se_adj = sqrt(sapply(1:k, 
                       adjusted_variance,
                       sigma2_vec = sigma2_vec, 
                       nu2 = nu2))
  
  CI_old = cbind("Unadjusted", 
                 theta_hat_vec + qnorm(1 - alpha/2)*se%*%t(c(-1, 1)))
  CrI_new = cbind("Adjusted", 
                  theta_hat_vec_adj + qnorm(1 - alpha/2)*se_adj%*%t(c(-1, 1)))
  CrI = as.data.frame(cbind(paste("Trt. Eff.", rep(1:k, 2)), 
                            rbind(CI_old, CrI_new)))
  names(CrI) = c("Effect", "Estimate", "LL", "UL")
  CrI = CrI %>% 
    mutate(LL = as.numeric(as.character(LL)),
           Center = c(theta_hat_vec, theta_hat_vec_adj),
           UL = as.numeric(as.character(UL)),
           Estimate = factor(Estimate, levels = c("Unadjusted", "Adjusted")),
           Effect = factor(Effect, levels = paste("Trt. Eff.", rep(1:k)))
    ) %>% 
    dplyr::filter(Estimate == 'Adjusted') %>% 
    mutate(`Est [95% CrI]` = paste0(sprintf('%.2f', Center), ' [',
                                    sprintf('%.2f', LL), ',',
                                    sprintf('%.2f', UL), ']')) %>% 
    select(`Est [95% CrI]`)
  
  return(CrI)
}



#*******************************************************************************
#* Inverting the covariance matrix of the marginal joint distribution of the 
#* treatment effect estimates
#* 
#* Inputs: 
#*   nu2 - the study-to-study variance
#*   sigma2_vec - the squared standard errors of the individual studies
#*   
#* Output: a matrix
#*******************************************************************************
Block_inverse = function(nu2, sigma2_vec){
  k = length(sigma2_vec)
  v = 1/(sigma2_vec + nu2)
  sum_v = sum(v)
  Mat = v%*%t(v)
  
  diag(v) - Mat/sum_v
}



#*******************************************************************************
#* The marginal log-likelihood of the data
#* 
#* Inputs: 
#*   log_nu2 - the study-to-study variance on the logarithmic scale
#*   sigma2_vec - the squared standard errors of the individual studies
#*   theta_hat_vec - a vector of treatment effect estimated before borrowing
#*   
#* Output: a real number
#*******************************************************************************
log_like = function(sigma2_vec, log_nu2, theta_hat_vec){
  nu2 = exp(log_nu2)
  v = 1/(sigma2_vec + nu2)
  Sigma_Inv = Block_inverse(nu2, sigma2_vec)
  
  - sum(log(v)) + log(sum(v)) + t(theta_hat_vec)%*%Sigma_Inv%*%theta_hat_vec
}



#*******************************************************************************
#* Empirical Bayes estimation of the study-to-study variance
#* 
#* Inputs: 
#*   sigma2_vec - the squared standard errors of the individual studies
#*   theta_hat_vec - a vector of treatment effect estimated before borrowing
#*   
#* Output: a real number
#*******************************************************************************
nu2_Emp_Bayes = function(sigma2_vec, theta_hat_vec){
  exp(nlminb(log(.5), log_like, 
             sigma2_vec = sigma2_vec,
             theta_hat_vec = theta_hat_vec)$par)
}




#*******************************************************************************
#* Study-specific Bayesian update of all treatment effects jointly, 
#* using empirical Bayes
#* 
#* Inputs: 
#*   theta_hat_vec - a vector of treatment effect estimated before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   alpha - 1 - the coverage level of the Bayesian credible intervals
#*   
#* Output: data frame containing treatment effect estimates and uncertainty 
#*         intervals, before and after Bayesian borrowing
#*******************************************************************************
output = function(theta_hat_vec, sigma2_vec, alpha){
  nu2_hat = nu2_Emp_Bayes(sigma2_vec, theta_hat_vec)
  
  k = length(theta_hat_vec)
  
  theta_hat_vec_adj = sapply(1:k, adjusted_estimate,
                             theta_hat_vec = theta_hat_vec, 
                             sigma2_vec = sigma2_vec, nu2 = nu2_hat)
  se = sqrt(sigma2_vec)
  
  se_adj = sqrt(sapply(1:k, 
                       adjusted_variance,
                       sigma2_vec = sigma2_vec, 
                       nu2 = nu2_hat))
  
  CI_old = cbind("Unadjusted", theta_hat_vec + qnorm(1 - alpha/2)*se%*%t(c(-1, 1)))
  CI_new = cbind("Adjusted", theta_hat_vec_adj + qnorm(1 - alpha/2)*se_adj%*%t(c(-1, 1)))
  CI = as.data.frame(cbind(paste("Trt.", rep(1:k, 2)), rbind(CI_old, CI_new)))
  names(CI) = c("Effect", "Estimate", "LL", "UL")
  CI$LL = as.numeric(as.character(CI$LL))
  CI$Center = c(theta_hat_vec, theta_hat_vec_adj)
  CI$UL = as.numeric(as.character(CI$UL))
  CI$Estimate = factor(CI$Estimate, levels = c("Unadjusted", "Adjusted"))
  CI$Effect = factor(CI$Effect, levels = paste("Trt.", rep(1:k)))
  
  return(list(CIs = CI, nu2_hat = nu2_hat))
}




#*******************************************************************************
#* Plotting the log-marginal likelihood vs. the study-to-study variance
#* 
#* Inputs: 
#*   sigma2_vec - the squared standard errors of the individual studies
#*   theta_hat_vec - a vector of treatment effect estimated before borrowing
#*   nu2_hat - the maximum likelihood estimate of nu2
#*   
#* Output: a ggplot
#*******************************************************************************
nu2_est_plot = function(sigma2_vec, theta_hat_vec, nu2_hat){
  nu2_vec = seq(0, nu2_hat*5, length = 100)
  loglikes = -sapply(log(nu2_vec), log_like, 
                     sigma2_vec = sigma2_vec, 
                     theta_hat_vec = theta_hat_vec)
  l_max = - log_like(sigma2_vec, log(nu2_hat), theta_hat_vec)
  likeli_df = data.frame(nu2 = nu2_vec, loglike = loglikes)
  
  ggplot(likeli_df, aes(nu2, loglike)) + 
    geom_line(linewidth = 1, col = 'dark red') + 
    geom_segment(aes(x = nu2_hat, xend = nu2_hat,
                     y = min(loglikes), 
                     yend = l_max), 
                 linetype = 'dashed', col = 'blue') + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    geom_point(aes(nu2_hat, l_max), 
               shape = 21, size = 3, fill = "white", 
               colour = "black", stroke = 2) + 
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=10),
          axis.line=element_line()) + 
    xlab(expression(nu^2)) + 
    ylab("Log-Marginal Likelihood")
}




#*******************************************************************************
#* A full Bayesian model run for the joint analysis of the treatment effects
#* 
#* Inputs: 
#*   theta_hat_vec - a vector of treatment effect estimated before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   N_iter - number of iterations per chain
#*   N_chains - number of chains
#*   thin_fact - the thinning factor
#*   
#* Output: a Stan object
#*******************************************************************************
stan_run = function(theta_hat_vec, sigma2_vec,
                    N_iter, N_chains, thin_fact){
  rstan_options(auto_write = TRUE)
  stan_comp = stan_model("Full_Bayes.stan")
  options(mc.cores = parallel::detectCores())
  fit = sampling(stan_comp, 
                 list(N = length(theta_hat_vec), 
                      theta_hat = theta_hat_vec, 
                      sigma = sqrt(sigma2_vec)),
                 iter = N_iter, chains = N_chains,
                 thin = thin_fact)
  
  return(fit)
}




#*******************************************************************************
#* Plotting a histogram of the posterior distribution of the study-to-study 
#* variance
#* 
#* Inputs: 
#*   fit - the output of the stan_run() function
#*   x_lim - a limit on the right tail of the plot
#*   
#* Output: a ggplot
#*******************************************************************************
nu2_posterior_plot = function(fit, x_lim){
  params = rstan::extract(fit)
  
  df = data.frame(nu2 = (params$nu)^2) %>% 
    filter(nu2 <= x_lim)
  
  d = density(df$nu2)
  ind = which.max(d$y)
  ymap = d$y[ind]
  xmap = d$x[ind]
  
  ggplot(df, aes(x = nu2)) +
    geom_histogram(df, mapping = aes(x = nu2, y = after_stat(density)),
                   fill = 'cadetblue2',
                   bins = 75, col = 'gray') +
    geom_density(size = 1, col = 'dark red', bw = "nrd0", 
                 kernel = "gaussian", na.rm = T) + 
    scale_x_continuous(expand = expansion(c(0,0.1)), limits = c(0, x_lim)) +
    scale_y_continuous(expand = c(0,0)) + 
    geom_segment(aes(x = xmap, xend = xmap, y = 0, yend = ymap),
                 linetype = 'dashed') + 
    geom_point(aes(xmap, ymap),
               shape = 21, size = 3, fill = 'white', stroke = 1) + 
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=10),
          axis.line=element_line()) +
    xlab(expression(nu^2)) 
}




#*******************************************************************************
#* Model diagnostics for the full Bayesian model
#* 
#* Inputs: 
#*   fit - the output of the stan_run() function
#*   
#* Output: a ggplot
#*******************************************************************************
stan_diagnostics = function(fit){
  p1 = stan_plot(fit, show_density = TRUE, ci_level = 0.95, 
                 fill_color = "cadetblue2", outer_level = .995, 
                 pars = 'theta') 
  
  p2 = stan_ac(fit, pars = c('theta', 'nu'), fill = 'cadetblue2') + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0))
  
  p3 = stan_trace(fit, pars = c('theta', 'nu'))
  
  (p1 + p2)/p3
}




#*******************************************************************************
#* Comparing the different nu2 estimation approaches with the raw estimate
#* 
#* Inputs: 
#*   stan_fit - the output of the stan_run() function
#*   out_Emp_Bayes - the output of the output() function
#*   theta_vec - the true treatment effects
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   out_fixed - the output of the output_fixed_nu() function
#*   
#* Output: a list consisting of a ggplot and a data frame
#*******************************************************************************
Method_comparison = function(stan_fit, out_Emp_Bayes,
                             theta_vec, theta_hat_vec,
                             out_fixed){
  CrI_MCMC = summary(stan_fit)$summary[1:length(theta_hat_vec), c(4,6,8)] %>%
    cbind(theta_vec) %>%
    apply(2, function(x){sprintf("%.2f", x)}) %>%
    as.data.frame %>%
    mutate(
      CrI = paste0(`50%`, " [", `2.5%`, ",", `97.5%`, "]")
    ) %>%
    dplyr::select(`Eff. size` = theta_vec, `Full Bayes [95% CrI]` = CrI)
  
  true_eff = CrI_MCMC[,1]
  CrI_MCMC = CrI_MCMC[,2]
  
  Emp_Bayes_CrI = out_Emp_Bayes$CIs %>%
    filter(Estimate == "Adjusted") %>%
    cbind(theta_vec) %>%
    select(LL, Center, UL, theta = theta_vec) %>%
    apply(2, as.character) %>%
    apply(2, as.numeric) %>%
    apply(2, function(x){sprintf("%.2f", x)}) %>%
    as.data.frame %>%
    mutate(
      CrI = paste0(Center, " [", LL, ",", UL, "]")
    ) %>%
    dplyr::select(`Emp. Bayes [95% CrI]` = CrI)
  
  
  Old_CrI = out_Emp_Bayes$CIs %>%
    filter(Estimate == "Unadjusted") %>%
    cbind(theta_vec) %>%
    dplyr::select(LL, Center, UL, theta = theta_vec) %>%
    apply(2, as.character) %>%
    apply(2, as.numeric) %>%
    apply(2, function(x){sprintf("%.2f", x)}) %>%
    as.data.frame %>%
    mutate(
      CrI = paste0(Center, " [", LL, ",", UL, "]")
    ) %>%
    dplyr::select(`Unadjusted [95% CI]` = CrI)
  
  
  tab = as.data.frame(cbind(true_eff, out_fixed, CrI_MCMC, 
                            Emp_Bayes_CrI, Old_CrI))
  names(tab) = c('Eff. size', 'Max BFP', 'Full Bayes', 
                 'Empirical Bayes', 'No Borrowing')
  rownames(tab) = c()
  
  helper = function(u){
    temp = sapply(unlist(u), function(x){strsplit(x, ' ')})
    Est = as.numeric(unlist(lapply(temp, `[[`, 1)))
    temp = lapply(lapply(temp, `[[`, 2), 
                  function(x){strsplit(x, ',')})
    temp = lapply(temp, `[[`, 1)
    LL = as.numeric(unlist(lapply(lapply(temp, `[[`, 1), 
                                  function(x){
                                    substr(x, 2, nchar(x))
                                  })))
    UL = as.numeric(unlist(lapply(lapply(temp, `[[`, 2), 
                                  function(x){
                                    substr(x, 1, nchar(x) - 1)
                                  })))
    
    return(data.frame(`Trt. Eff.` = paste0('$\\theta_', 
                                           1:length(Est),
                                           '$'), 
                      Est = Est, LL = LL, UL = UL,
                      check.names = F))
  }
  
  df_Naive = cbind(Model = 'No Borrowing', helper(Old_CrI))
  df_Emp = cbind(Model = 'Empirical Bayes', helper(Emp_Bayes_CrI))
  df_MCMC = cbind(Model = 'Full Bayes', helper(CrI_MCMC))
  df_Fixed = cbind(Model = 'Max BFP', helper(out_fixed))
  
  df = rbind(df_Naive, df_Emp, df_MCMC, df_Fixed) %>% 
    mutate(True = as.numeric(rep(true_eff, 4)))
  
  p = ggplot(df, aes(`Trt. Eff.`, Est, col = Model, 
                     bg = Model, shape = Model)) + 
    geom_errorbar(df, mapping = aes(ymin = LL, ymax = UL),
                  position = position_dodge(width = 0.5), width = .2) + 
    geom_point(position = position_dodge(width = 0.5), size = 3) + 
    scale_x_discrete(labels = latex2exp::TeX(unique(df$`Trt. Eff.`))) + 
    coord_flip() +
    theme(panel.background = element_blank(),
          legend.position = 'top',
          panel.grid.major = element_line(colour = "grey92"),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=18, face="bold"),
          axis.title.y = element_text(size=18, face="bold"),
          axis.text.y = element_text(size=14),
          legend.text=element_text(size=12), 
          axis.line=element_line()) + 
    scale_fill_brewer(palette = "Set1") + 
    scale_color_brewer(palette = "Set1") + 
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    geom_segment(df, 
                 mapping = aes(x = `Trt. Eff.`, 
                               xend = `Trt. Eff.`,
                               y = True, yend = True),
                 position = position_dodge(width = 0.75),
                 linetype = 'dashed') + 
    ylab('Estimate') + 
    scale_shape_manual(values = c(21:24))
  
  return(list(Table = tab, Plot = p))
}




#*******************************************************************************
#* Calculating the theoretical mean squared errors
#* 
#* Inputs: 
#*   nu2 - the study-to-study variance
#*   sigma2_vec - the squared standard errors of the individual studies
#*   theta_vec - the true treatment effects
#*   i - the index of the target study to be updated by borrowing information 
#*       from other studies
#*   
#* Output: a real number
#*******************************************************************************
MSE = function(nu2, sigma2_vec, theta_vec, i){
  K = length(sigma2_vec)
  
  theta = theta_vec[-i]
  theta_diff = theta - theta_vec[i]
  
  s2_k = sigma2_vec[i]
  u = 1/(sigma2_vec + nu2)
  
  denom = (sum(u))^2
  numer = (t(u[-i])%*%theta_diff)^2*s2_k + (1 + nu2*sum(u[-i]))^2
  
  s2_k/(s2_k + nu2)^2*numer/denom
}



#*******************************************************************************
#* A wrapper of MSE() that receives (nu2, i) as a vector
#*******************************************************************************
MSE_vec = function(vec, sigma2_vec, theta_vec){
  MSE(vec[1], sigma2_vec, theta_vec, vec[2])
}




#*******************************************************************************
#* Plotting the theoretical root mean squared errors vs. the study-tostudy 
#* variance for all treatment effects
#* 
#* Inputs: 
#*   nu2_vec - a vector of nu2 values for RMSE evaluation
#*   sigma2_vec - the squared standard errors of the individual studies
#*   nu2_hat - the MLE of the study-to-study variance (for highlighting)
#*   theta_vec - the true treatment effects
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   
#* Output: a list consisting of a ggplot and a data frame
#*******************************************************************************
RMSE_vs_nu2_plot_table = function(nu2_vec, sigma2_vec, nu2_hat,
                                  theta_vec, theta_hat_vec){
  df_adjust = data.frame(nu2 = rep(nu2_vec, length(sigma2_vec)), 
                         param = rep(1:length(sigma2_vec),
                                     each = length(nu2_vec)),
                         Model = 'Informed')
  
  df_adjust$RMSE = sqrt(apply(df_adjust[,1:2], 1, MSE_vec,
                              sigma2_vec = sigma2_vec, 
                              theta_vec = theta_vec)) 
  
  df_adjust_avg = df_adjust %>% 
    group_by(nu2) %>% 
    summarise(RMSE = sqrt(mean(RMSE^2))) %>% 
    mutate(param = 'Average',
           Model = 'Informed') %>% 
    relocate(RMSE, .after = Model) %>% 
    as.data.frame()
  
  df_adjust = rbind(df_adjust, df_adjust_avg)
  
  df_unadjust = data.frame(nu2 = rep(nu2_vec, length(sigma2_vec)),
                           param = rep(1:length(sigma2_vec),
                                       each = length(nu2_vec)),
                           Model = 'Uninformed',
                           RMSE = sqrt(rep(sigma2_vec, 
                                           each = length(nu2_vec))))
  
  df_unadjust_avg = df_unadjust %>% 
    group_by(nu2) %>% 
    summarise(RMSE = sqrt(mean(RMSE^2))) %>% 
    mutate(param = 'Average',
           Model = 'Uninformed') %>% 
    relocate(RMSE, .after = Model) %>% 
    as.data.frame()
  
  df_unadjust = rbind(df_unadjust, df_unadjust_avg)
  
  df = rbind(df_adjust, df_unadjust) %>% 
    mutate(param = ifelse(param != 'Average', 
                          paste0('theta[', param, ']'),
                          param)) %>% 
    mutate(param = factor(param, 
                          levels = c(paste0('theta[',
                                            1:length(sigma2_vec),
                                            ']'),
                                     'Average'))
    )
  
  out = output(theta_hat_vec, sigma2_vec, alpha = .05)
  nu2_hat = out$nu2_hat 
  
  df_line = df %>% 
    filter(abs(nu2 - nu2_hat) == min(abs(nu2 - nu2_hat))) %>% 
    select(-nu2) %>% 
    reshape2::dcast(param ~ Model) %>% 
    mutate(Model = 'Informed')
  
  
  p = ggplot(df, aes(nu2, RMSE, col = Model, linetype = Model)) + 
    geom_line(size = .75) + 
    scale_y_continuous(expand = expansion(c(0, 0.1))) + 
    scale_x_continuous(expand = expansion(c(0, 0.025))) + 
    facet_rep_wrap(.~ param,
                   labeller = label_parsed,
                   repeat.tick.labels = 'all') + 
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          legend.position = "top",
          legend.text = element_text(size = 12), 
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(hjust = 1, size = 12),
          plot.title = element_text(size = 16, face = "bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.line=element_line()) + 
    scale_color_brewer(palette="Set1") + 
    xlab(expression(nu^2)) + 
    geom_segment(df_line, 
                 mapping = aes(x = nu2_hat, xend = nu2_hat,
                               y = Uninformed, yend = Informed),
                 linetype = 'dashed', col = 1) + 
    geom_point(df_line,
               mapping = aes(x = nu2_hat, y = Informed), 
               shape = 21, size = 1.75, fill = "white", 
               colour = "black", stroke = 1)
  
  df_out = df_line %>% 
    select(-Model) %>% 
    mutate(Informed = sprintf('%.3f', Informed),
           Uninformed = sprintf('%.3f', Uninformed)) %>% 
    relocate(Informed, .after = Uninformed) %>% 
    rename('Parameter' = 'param',
           'Informed RMSE' = 'Informed',
           'Uninformed RMSE' = 'Uninformed') %>% 
    mutate(Parameter = paste0('$\\theta_', 1:nrow(df_line), 
                              ' = ', sprintf('%.2f', theta_vec), '$')) %>% 
    kable(align = 'c') %>%
    kable_styling(full_width = F, font_size = 14,
                  bootstrap_options = c("striped")) %>%
    row_spec(0, font_size = 14) 
  
  
  return(list(Plot = p, Table = df_out))
}



#*******************************************************************************
#* Standard error for the log-relative risk estimator
#* 
#* Inputs: 
#*   n - number of patients per arm
#*   p - the control event rate
#*   theta - the (true) log-relative risk
#*   
#* Output: a number
#*******************************************************************************
trt_eff_to_sd = function(n, p, theta){
  sqrt(1/n*(1/p*(1 + exp(-theta)) - 2))
}



#*******************************************************************************
#* Standard error for the log-odds ratio estimator
#* 
#* Inputs: 
#*   n - number of patients per arm
#*   p - the control event rate
#*   theta - the (true) log-odds ratio
#*   
#* Output: a number
#*******************************************************************************
trt_eff_to_sd_OR = function(n, p, theta){
  u = p*exp(theta)
  sig2 = 1/n*((2-p) + 1/u*(1 - p + u)^2)/(1 - p)
  
  sqrt(sig2)
}




#*******************************************************************************
#* Calculating the power of the hypothesis test based on the posterior 
#* probability of efficacy (using empirical Bayes to estimate nu2) for relative 
#* risks
#* 
#* Inputs: 
#*   n - number of patients per arm
#*   CER - the control event rate
#*   RR - the relative risk
#*   N_MC - number of Monte Carlo draws
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   alpha - 1 - the posterior probability of efficacy threshold
#*   
#* Output: a number
#*******************************************************************************
power_func_sim = function(n, CER, RR, N_MC, 
                          theta_hat_vec, sigma2_vec, alpha){
  set.seed(n*N_MC)
  
  y_ctrl = rbinom(N_MC, n, CER)
  y_trt = rbinom(N_MC, n, CER*RR)
  p_ctrl_hat = y_ctrl/n
  p_trt_hat = y_trt/n
  
  theta_hat_new = log(p_trt_hat/p_ctrl_hat)
  sigma2_hat_new = (unlist(trt_eff_to_sd(n, p_ctrl_hat, theta_hat_new)))^2
  
  df = cbind(theta_hat_new, sigma2_hat_new)
  
  power_aux_func = function(v){
    theta_hat_last = v[1]
    sigma2_hat_last = v[2]
    
    K = length(theta_hat_vec)
    
    sigma2_vec_aug = c(sigma2_vec, sigma2_hat_last)
    theta_hat_vec_aug = c(theta_hat_vec, theta_hat_last)
    
    nu2_hat = nu2_Emp_Bayes(sigma2_vec_aug, theta_hat_vec_aug)
    theta_hat_star = adjusted_estimate(theta_hat_vec_aug, sigma2_vec_aug, 
                                       nu2_hat, K+1)
    sigma_star = sqrt(adjusted_variance(sigma2_vec_aug, nu2_hat, K+1))
    
    post_effic = pnorm(-theta_hat_star/sigma_star)
    return(as.numeric(post_effic > (1 - alpha)))
  }
  
  
  pow = apply(df, 1, power_aux_func)
  return(mean(pow))
}



#*******************************************************************************
#* Calculating the power of the hypothesis test based on the posterior 
#* probability of efficacy (using empirical Bayes to estimate nu2) for odds 
#* ratios
#* 
#* Inputs: 
#*   n - number of patients per arm
#*   CER - the control event rate
#*   OR - the odds ratio
#*   N_MC - number of Monte Carlo draws
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   alpha - 1 - the posterior probability of efficacy threshold
#*   
#* Output: a number
#*******************************************************************************
power_func_sim_OR = function(n, CER, OR, N_MC, 
                             theta_hat_vec, sigma2_vec, alpha){
  set.seed(n*N_MC)
  
  p_trt = CER*OR/(1 - CER + CER*OR)
  
  y_ctrl = rbinom(N_MC, n, CER)
  y_trt = rbinom(N_MC, n, p_trt)
  p_ctrl_hat = y_ctrl/n
  p_trt_hat = y_trt/n
  
  theta_hat_new = log(p_trt_hat*(1 - p_ctrl_hat)/(p_ctrl_hat*(1 - p_trt_hat)))
  sigma2_hat_new = (unlist(trt_eff_to_sd_OR(n, p_ctrl_hat, theta_hat_new)))^2
  
  df = cbind(theta_hat_new, sigma2_hat_new)
  
  power_aux_func = function(v){
    theta_hat_last = v[1]
    sigma2_hat_last = v[2]
    
    K = length(theta_hat_vec)
    
    sigma2_vec_aug = c(sigma2_vec, sigma2_hat_last)
    theta_hat_vec_aug = c(theta_hat_vec, theta_hat_last)
    
    nu2_hat = nu2_Emp_Bayes(sigma2_vec_aug, theta_hat_vec_aug)
    theta_hat_star = adjusted_estimate(theta_hat_vec_aug, sigma2_vec_aug, 
                                       nu2_hat, K+1)
    sigma_star = sqrt(adjusted_variance(sigma2_vec_aug, nu2_hat, K+1))
    
    post_effic = pnorm(theta_hat_star/sigma_star)
    return(as.numeric(post_effic > (1 - alpha)))
  }
  
  
  pow = apply(df, 1, power_aux_func)
  return(mean(pow))
}



#*******************************************************************************
#* Calculating the power of the hypothesis test based on the posterior 
#* probability of efficacy (assuming a known nu2) for odds ratios
#* 
#* Inputs: 
#*   n - number of patients per arm
#*   CER - the control event rate
#*   OR - the odds ratio
#*   N_MC - number of Monte Carlo draws
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   alpha - 1 - the posterior probability of efficacy threshold
#*   nu2 - the study-to-study variance
#*   
#* Output: a number
#*******************************************************************************
power_func_sim_OR_fixed_nu2 = function(n, CER, OR, N_MC, 
                                       theta_hat_vec, sigma2_vec, 
                                       alpha, nu2){
  set.seed(n*N_MC)
  
  p_trt = CER*OR/(1 - CER + CER*OR)
  
  y_ctrl = rbinom(N_MC, n, CER)
  y_trt = rbinom(N_MC, n, p_trt)
  p_ctrl_hat = y_ctrl/n
  p_trt_hat = y_trt/n
  
  theta_hat_new = log(p_trt_hat*(1 - p_ctrl_hat)/(p_ctrl_hat*(1 - p_trt_hat)))
  sigma2_hat_new = (unlist(trt_eff_to_sd_OR(n, p_ctrl_hat, theta_hat_new)))^2
  
  df = cbind(theta_hat_new, sigma2_hat_new)
  
  power_aux_func = function(v){
    theta_hat_last = v[1]
    sigma2_hat_last = v[2]
    
    K = length(theta_hat_vec)
    
    sigma2_vec_aug = c(sigma2_vec, sigma2_hat_last)
    theta_hat_vec_aug = c(theta_hat_vec, theta_hat_last)
    
    theta_hat_star = adjusted_estimate(theta_hat_vec_aug, sigma2_vec_aug, 
                                       nu2, K+1)
    sigma_star = sqrt(adjusted_variance(sigma2_vec_aug, nu2, K+1))
    
    post_effic = pnorm(theta_hat_star/sigma_star)
    return(as.numeric(post_effic > (1 - alpha)))
  }
  
  
  pow = apply(df, 1, power_aux_func)
  return(mean(pow))
}



#*******************************************************************************
#* Calculating the posterior probability of efficacy in a beta-binomial model
#* 
#* Inputs: 
#*   n_events - a vector of two numbers of events (control, treatment)
#*   n_patients - a vector of two numbers of patients (control, treatment)
#*   N_MC - number of Monte Carlo draws
#*   
#* Output: a number
#*******************************************************************************
beta_post_efficacy = function(n_events, n_patients, N_MC){
  alpha = 1 + n_events
  beta = 1 + n_patients - n_events
  
  samp = matrix(rbeta(N_MC*2, 
                      rep(alpha, each = N_MC), 
                      rep(beta, each = N_MC)),
                ncol = 2)
  
  post_effic = mean(samp[,2] < samp[,1])
  
  return(post_effic)
}



#*******************************************************************************
#* A single beta-binomial simulation
#* 
#* Inputs: 
#*   CER - the control event rate
#*   RR - the relative risk
#*   n - number of patients per arm
#*   N_MC - number of Monte Carlo draws from the beta posterior distribution
#*   alpha - 1 - the posterior probability of efficacy threshold
#*   sim_num - the simulation seed
#*   
#* Output: a logical (null rejected or not)
#*******************************************************************************
beta_single_simulation = function(CER, RR, n, N_MC, alpha, sim_num){
  set.seed(sim_num)
  
  n_events = rbinom(2, n, CER*c(1, RR))
  post_eff = beta_post_efficacy(n_events, rep(n, 2), N_MC)
  
  return(post_eff > 1 - alpha)
}



#*******************************************************************************
#* Calculating the power of the hypothesis test based on the posterior 
#* probability of efficacy in a beta-binomial model
#* 
#* Inputs: 
#*   CER - the control event rate
#*   R - the relative risk
#*   n - number of patients per arm
#*   N_MC - number of Monte Carlo draws from the beta posterior distribution
#*   alpha - 1 - the posterior probability of efficacy threshold
#*   N_sim - number of random data sets to be simulated
#*   
#* Output: a number
#*******************************************************************************
beta_power_func = function(CER, RR, n, N_MC, alpha, N_sim){
  out = sapply(1:N_sim, 
               beta_single_simulation,
               CER = CER, RR = RR, n = n, 
               N_MC = N_MC, alpha = alpha)
  
  mean(out)
}



#*******************************************************************************
#* Calculating the power of the borrowing model against that of the 
#* beta-binomial model vs. the sample size (per arm)
#* 
#* Inputs: 
#*   n_vec - a vector of numbers of patients per arm for power evaluation
#*   CER - the control event rate
#*   N_MC - number of Monte Carlo draws from the beta posterior distribution
#*          (for the beta-binomial model)
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   alpha - 1 - the posterior probability of efficacy threshold
#*   N_sim - number of random data sets to be simulated
#*   
#* Output: a number
#*******************************************************************************
power_vs_n_df = function(n_vec, CER, RR, N_MC, theta_hat_vec, 
                         sigma2_vec, alpha, N_sim){
  cl = makeSOCKcluster(parallel::detectCores() - 1)
  registerDoSNOW(cl)
  
  pb = txtProgressBar(max = length(n_vec), style=3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  
  l = foreach(i = 1:length(n_vec),
              .options.snow = opts,
              .export = c(
                'nu2_Emp_Bayes',
                "trt_eff_to_sd",
                'log_like',
                'adjusted_estimate',
                'adjusted_variance',
                'a_coeffs',
                'power_func_sim',
                'Block_inverse',
                'beta_post_efficacy',
                'beta_single_simulation',
                'beta_power_func'),
              .packages = c('randtoolbox')) %dopar% {
                set.seed(i)
                n = n_vec[i]
                
                out = cbind(power_func_sim(n = n, CER = CER, RR = RR, 
                                           N_MC = N_sim, 
                                           theta_hat_vec = theta_hat_vec, 
                                           sigma2_vec = sigma2_vec, 
                                           alpha = alpha),
                            beta_power_func(CER, RR, n, N_MC, 
                                            alpha, N_sim))
                
              }
  close(pb)
  stopCluster(cl)
  registerDoSEQ()
  
  df = as.data.frame(cbind(n_vec, do.call(rbind, l)))
  names(df) = c('n', 'Informed', 'Uninformed')
  
  df = df %>% 
    melt(id.vars = 'n',
         variable.name = "Model",
         value.name = "Power")
  
  return(df)
}



#*******************************************************************************
#* Plotting the power of the borrowing model against that of the 
#* beta-binomial model vs. the sample size (per arm)
#* 
#* Inputs: 
#*   df - the output of the power_vs_n_df() function
#*   power_thresh - what power do we wish to highlight (e.g., 80%)
#*   
#* Output: a ggplot
#*******************************************************************************
power_vs_n_plot = function(df, power_thresh){
  df1 = df %>% filter(Model == 'Uninformed')
  
  fit1 = cobs(df1$n, df1$Power,
              constraint= "increase", 
              lambda=0, 
              degree=1, # for L1 roughness
              knots=seq(min(df1$n), max(df1$n),length.out=20), # desired nr of knots 
              tau=0.5) 
  df1$Power = predict(fit1, interval = "none", z = df1$n)[,2]
  
  df2 = df %>% filter(Model == 'Informed')
  
  fit2 = cobs(df2$n, df2$Power,
              constraint= "increase", 
              lambda=0, 
              degree=1, # for L1 roughness
              knots=seq(min(df2$n), max(df2$n),length.out=20), # desired nr of knots 
              tau=0.5) 
  df2$Power = predict(fit2, interval = "none", z = df2$n)[,2]
  
  df = rbind(df1, df2)
  
  df_line = df %>% 
    filter(Power >= power_thresh) %>% 
    group_by(Model)  %>% 
    filter(row_number() == 1)
  
  print(df_line %>% 
          mutate(Power = paste0(sprintf('%.1f', Power*100), '%')) %>% 
          as.data.frame)
  
  
  ggplot(df, aes(n, Power, col = Model)) + 
    geom_line(linewidth = 1)  +  
    scale_x_continuous(expand = expansion(mult = c(0, 0.025))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.025)),
                       limits = c(0, 1), labels = percent, 
                       breaks = .2*c(0:5)) + 
    geom_segment(df_line, 
                 mapping = aes(x = min(df$n), xend = n,
                               y = power_thresh, yend = power_thresh),
                 linetype = 'dashed') + 
    geom_segment(df_line, 
                 mapping = aes(x = n, xend = n,
                               y = 0, yend = power_thresh),
                 linetype = 'dashed') + 
    geom_point(df_line, mapping = aes(n, Power), shape = 21, 
               size = 2, fill = 'white', stroke = 2) + 
    theme(panel.background = element_blank(),
          legend.position = "top",
          panel.grid.major = element_line(colour = "grey92"),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=10),
          axis.line=element_line()) + 
    scale_color_brewer(palette="Set1") + 
    xlab('Sample size (per arm)')
}




#*******************************************************************************
#* Calculating the power of the borrowing model against that of the 
#* beta-binomial model vs. the relative risk reduction (for a fixed sample size)
#* 
#* Inputs: 
#*   n - the sample size (per arm)
#*   CER - the control event rate
#*   RR_vec - a vector of relative risks for power evaluation
#*   N_MC - number of Monte Carlo draws from the beta posterior distribution
#*          (for the beta-binomial model)
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   alpha - 1 - the posterior probability of efficacy threshold
#*   N_sim - number of random data sets to be simulated
#*   
#* Output: a number
#*******************************************************************************
power_vs_RRR_df = function(n, CER, RR_vec, N_MC, theta_hat_vec, 
                           sigma2_vec, alpha, N_sim){
  cl = makeSOCKcluster(parallel::detectCores() - 1)
  registerDoSNOW(cl)
  
  pb = txtProgressBar(max = length(RR_vec), style=3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress=progress)
  
  l = foreach(i = 1:length(RR_vec),
              .options.snow = opts,
              .export = c(
                'nu2_Emp_Bayes',
                "trt_eff_to_sd",
                'log_like',
                'adjusted_variance',
                'a_coeffs',
                'adjusted_estimate',
                'power_func_sim',
                'Block_inverse',
                'beta_post_efficacy',
                'beta_single_simulation',
                'beta_power_func'),
              .packages = c('randtoolbox')) %dopar% {
                set.seed(i)
                RR = RR_vec[i]
                
                out = cbind(power_func_sim(n = n, CER = CER, RR = RR, 
                                           N_MC = N_sim, 
                                           theta_hat_vec = theta_hat_vec, 
                                           sigma2_vec = sigma2_vec, 
                                           alpha = alpha),
                            beta_power_func(CER, RR, n, N_MC, 
                                            alpha, N_sim))
                
              }
  close(pb)
  stopCluster(cl)
  registerDoSEQ()
  
  df = as.data.frame(cbind(1 - RR_vec, do.call(rbind, l)))
  names(df) = c('RRR', 'Informed', 'Uninformed')
  
  df = df %>% 
    melt(id.vars = 'RRR',
         variable.name = "Model",
         value.name = "Power")
  
  return(df)
}




#*******************************************************************************
#* Calculating the power of the borrowing model against that of the 
#* beta-binomial model vs. the odds ratio (for a fixed sample size)
#* 
#* Inputs: 
#*   n - the sample size (per arm)
#*   CER - the control event rate
#*   OR_vec - a vector of odds ratios for power evaluation
#*   N_MC - number of Monte Carlo draws from the beta posterior distribution
#*          (for the beta-binomial model)
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   alpha - 1 - the posterior probability of efficacy threshold
#*   N_sim - number of random data sets to be simulated
#*   
#* Output: a number
#*******************************************************************************
power_vs_OR_df = function(n, CER, OR_vec, N_MC, theta_hat_vec, 
                          sigma2_vec, alpha, N_sim){
  cl = makeSOCKcluster(parallel::detectCores() - 1)
  registerDoSNOW(cl)
  
  pb = txtProgressBar(max = length(OR_vec), style=3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress=progress)
  
  l = foreach(i = 1:length(RR_vec),
              .options.snow = opts,
              .export = c(
                'nu2_Emp_Bayes',
                "trt_eff_to_sd",
                'log_like',
                'adjusted_variance',
                'a_coeffs',
                'adjusted_estimate',
                'power_func_sim',
                'Block_inverse',
                'beta_post_efficacy',
                'beta_single_simulation',
                'beta_power_func'),
              .packages = c('randtoolbox')) %dopar% {
                set.seed(i)
                OR = OR_vec[i]
                
                out = power_func_sim_OR(n = n, CER = CER, OR = OR, 
                                        N_MC = N_sim, 
                                        theta_hat_vec = theta_hat_vec, 
                                        sigma2_vec = sigma2_vec, 
                                        alpha = alpha)
                
              }
  close(pb)
  stopCluster(cl)
  registerDoSEQ()
  
  df = as.data.frame(cbind(OR_vec, do.call(rbind, l)))
  names(df) = c('OR', 'Power')
  
  return(df)
}



#*******************************************************************************
#* Plotting the power of the borrowing model against that of the 
#* beta-binomial model vs. the relative risk reduction (for a fixed sample size)
#* 
#* Inputs: 
#*   df - the output of the power_vs_RRR_df() function
#*   power_thresh - what power do we wish to highlight (e.g., 80%)
#*   
#* Output: a number
#*******************************************************************************
power_vs_RRR_plot = function(df, power_thresh){
  df1 = df %>% filter(Model == 'Uninformed')
  
  fit1 = cobs(df1$RRR, df1$Power,
              constraint= "increase", 
              lambda=0, 
              degree=1, # for L1 roughness
              knots=seq(min(df1$RRR), max(df1$RRR),length.out=50), # desired nr of knots 
              tau=0.5) 
  df1$Power = predict(fit1, interval = "none", z = df1$RRR)[,2]
  
  df2 = df %>% filter(Model == 'Informed')
  
  fit2 = cobs(df2$RRR, df2$Power,
              constraint= "increase", 
              lambda=0, 
              degree=1, # for L1 roughness
              knots=seq(min(df2$RRR), max(df2$RRR),length.out=50), # desired nr of knots 
              tau=0.5) 
  df2$Power = predict(fit2, interval = "none", z = df2$RRR)[,2]
  
  df = rbind(df1, df2)
  
  df_line = df %>% 
    filter(Power >= power_thresh) %>% 
    group_by(Model)  %>% 
    arrange(RRR) %>% 
    filter(row_number() == 1)
  
  print(df_line %>% 
          mutate(RRR = paste0(sprintf('%.1f', RRR*100), '%'),
                 Power = paste0(sprintf('%.1f', Power*100), '%')) %>% 
          as.data.frame)
  
  
  ggplot(df, aes(RRR, Power, col = Model)) + 
    geom_line(linewidth = 1)  +  
    scale_x_continuous(expand = expansion(mult = c(0, 0.025)),
                       labels = percent) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.025)),
                       limits = c(0, 1), labels = percent, 
                       breaks = .2*c(0:5)) + 
    geom_segment(df_line, 
                 mapping = aes(x = min(df$RRR), xend = RRR,
                               y = power_thresh, yend = power_thresh),
                 linetype = 'dashed') + 
    geom_segment(df_line, 
                 mapping = aes(x = RRR, xend = RRR,
                               y = 0, yend = power_thresh),
                 linetype = 'dashed') + 
    geom_point(df_line, mapping = aes(RRR, Power), shape = 21, 
               size = 2, fill = 'white', stroke = 2) + 
    theme(panel.background = element_blank(),
          legend.position = "top",
          panel.grid.major = element_line(colour = "grey92"),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=10),
          axis.line=element_line()) + 
    scale_color_brewer(palette="Set1") + 
    xlab('Relative risk reduction')
}





#*******************************************************************************
#* A single simulation for the assessment of Bayesian borrowing performance vs.
#* no borrowing (using empirical Bayes to estimate nu2)
#* 
#* Inputs: 
#*   K - the number of studies
#*   p_ctrl - control event rate
#*   p_trt - treatment event rate
#*   sigma2_ctrl - variance for the control group study-to-study perturbation
#*   sigma2_trt - variance for the treatment group study-to-study perturbation
#*   n_ctrl - control group number of patients (per study)
#*   n_trt - treatment group number of patients (per study)
#*   alpha - 1 - uncertainty intervals coverage
#*   sim_num - simulation seed
#*   
#* Output: a two-row data frame consisting of estimation accuracy and coverage 
#*         assessment for the informed and uninformed models in a single random 
#*         dataset
#*******************************************************************************
single_simulation = function(K, p_ctrl, p_trt, sigma2_ctrl, sigma2_trt, 
                             n_ctrl, n_trt, alpha, sim_num){
  set.seed(sim_num)
  
  dat = random_data(K = K, 
                    p_ctrl = p_ctrl, 
                    p_trt = p_trt,
                    sigma2_ctrl = sigma2_ctrl, 
                    sigma2_trt = sigma2_trt,
                    n_ctrl = n_ctrl, n_trt = n_trt)
  
  out = output(dat$theta_hat_vec, dat$sigma2_vec, alpha = alpha)
  CI = out$CIs
  
  
  unadj = cbind(CI[1:K,c(1,3:5)], dat$theta_vec)
  adj = cbind(CI[(K+1):(2*K),c(1,3:5)], dat$theta_vec)
  
  bias_SqE_cover = function(subdat){
    delta = subdat[,5] - subdat[,4]
    SqE = delta^2
    CI_length = subdat[,3] - subdat[,2]
    covered = (subdat[,5] - subdat[,2])*(subdat[,5] - subdat[,3]) <= 0
    return(cbind(sim_num, Effect = subdat$Effect, delta, 
                 SqE, CI_length, covered))
  }
  
  out = rbind(cbind(Type = 'Unadjusted', 
                    as.data.frame(bias_SqE_cover(unadj))),
              cbind(Type = 'Adjusted', 
                    as.data.frame(bias_SqE_cover(adj)))
  )
  
  return(out)
}




#*******************************************************************************
#* A full simulation for the assessment of Bayesian borrowing performance vs.
#* no borrowing (using empirical Bayes to estimate nu2)
#* 
#* Inputs: 
#*   K - the number of studies
#*   p_ctrl - control event rate
#*   p_trt - treatment event rate
#*   sigma2_ctrl - variance for the control group study-to-study perturbation
#*   sigma2_trt - variance for the treatment group study-to-study perturbation
#*   n_ctrl - control group number of patients (per study)
#*   n_trt - treatment group number of patients (per study)
#*   alpha - 1 - uncertainty intervals coverage
#*   N_sims - number of simulations
#*   
#* Output: a large data frame consisting of estimation accuracy and coverage 
#*         assessment for the informed and uninformed models over a large number
#*         of simulations
#*******************************************************************************
adj_vs_unadj_simulation = function(K, p_ctrl, p_trt, sigma2_ctrl, sigma2_trt, 
                                   n_ctrl, n_trt, alpha, N_sims){
  cl = parallel::makeCluster(detectCores() - 1)
  registerDoSNOW(cl)
  
  pb = txtProgressBar(max = N_sims, style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  
  sout = foreach(s = 1:N_sims,
                 .options.snow = opts,
                 .export = c(
                   "random_data",
                   'nu2_Emp_Bayes',
                   "output",
                   'single_simulation',
                   'log_like',
                   'sigma2_vec',
                   'Block_inverse',
                   'theta_hat_vec',
                   'adjusted_estimate',
                   'adjusted_variance'
                 ),
                 .packages = c('LaplacesDemon')) %dopar% {
                   out = single_simulation(K, p_ctrl, p_trt, sigma2_ctrl, 
                                           sigma2_trt, n_ctrl, n_trt, alpha, 
                                           sim_num = s)
                   
                   return(out)
                 }
  
  close(pb)
  stopCluster(cl)
  registerDoSEQ()
  
  return(do.call(rbind, sout))
}



#*******************************************************************************
#* Converting a Bayesian borrowing fraction into a study-to-study variance
#* 
#* Inputs: 
#*   sigma2_vec - the squared standard errors of the individual studies
#*   target_weight - the Bayesian borrowing fraction to be converted
#*   
#* Output: a number
#*******************************************************************************
weight_to_nu2 = function(sigma2_vec, target_weight){
  
  weights_diff = function(nu2){
    weights(sigma2_vec, nu2, 3)[1:2] - target_weight
  }
  
  nu2_seq = seq(1e-6, .1, by = 1e-6)
  temp = t(sapply(nu2_seq, weights_diff))
  
  nu2_target = nu2_seq[max(apply(abs(temp), 2, which.min))]
  return(nu2_target)
}




#*******************************************************************************
#* Bayesian posterior inference with a known borrowing fraction
#* 
#* Inputs: 
#*   theta_hat_vec - the estimated treatment effects before borrowing
#*   sigma2_vec - the squared standard errors of the individual studies
#*   alpha - 1 - the posterior probability of efficacy threshold
#*   w - the borrowing fraction used
#*   
#* Output: posterior mean and credible interval
#*******************************************************************************
output_fixed_w = function(theta_hat_vec, sigma2_vec, alpha, w){
  nu2 = weight_to_nu2(sigma2_vec, w)
  
  k = length(theta_hat_vec)
  theta_hat_vec_adj = c(adjusted_estimate(theta_hat_vec, 
                                          sigma2_vec, 
                                          nu2, k))
  
  se_adj = sqrt(adjusted_variance(sigma2_vec, nu2, k))
  
  c(theta_hat_vec_adj, theta_hat_vec_adj + 
      c(-1,1)*qnorm(1-alpha/2)*se_adj)
}