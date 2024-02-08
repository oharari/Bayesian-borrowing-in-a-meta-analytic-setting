random_data = function(K, p_ctrl, p_trt, sigma2_ctrl, 
                       sigma2_trt, n_ctrl, n_trt){
  p_ctrl_vec = invlogit(rnorm(K, logit(p_ctrl), sigma2_ctrl))
  p_trt_vec = invlogit(rnorm(K, logit(p_trt), sigma2_trt))
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




weights = function(sigma2_vec, nu2, i){
  v = 1/(sigma2_vec + nu2)
  u = sigma2_vec[i]*v
  w = sigma2_vec/sigma2_vec[i]*u/(u[i] + nu2*sum(v))
  w[i] = 1
  
  w
}




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





adjusted_estimate = function(theta_hat_vec, 
                             sigma2_vec, nu2, i){
  k = length(sigma2_vec)
  v = 1/(sigma2_vec + nu2)
  sum_v = sum(v)
  a = sigma2_vec[i]*v[i]*v
  a[i] = a[i] + nu2*v[i]*sum_v
  
  return(a%*%theta_hat_vec/sum_v)
}





adjusted_variance = function(sigma2_vec, nu2, i){
  v = 1/(sigma2_vec + nu2)
  sum_v = sum(v)
  term = sigma2_vec[i]*v[i]
  
  term*(term + nu2*sum_v)/sum_v
}





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




Block_inverse = function(nu2, sigma2_vec){
  k = length(sigma2_vec)
  v = 1/(sigma2_vec + nu2)
  sum_v = sum(v)
  Mat = v%*%t(v)
  
  diag(v) - Mat/sum_v
}




log_like = function(sigma2_vec, log_nu2, theta_hat_vec){
  nu2 = exp(log_nu2)
  v = 1/(sigma2_vec + nu2)
  Sigma_Inv = Block_inverse(nu2, sigma2_vec)
  
  - sum(log(v)) + log(sum(v)) + t(theta_hat_vec)%*%Sigma_Inv%*%theta_hat_vec
}




nu2_Emp_Bayes = function(sigma2_vec, theta_hat_vec){
  exp(nlminb(log(.5), log_like, 
             sigma2_vec = sigma2_vec,
             theta_hat_vec = theta_hat_vec)$par)
}




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




nu2_est_plot = function(sigma2_vec, theta_hat_vec, nu2_hat){
  nu2_vec = seq(0, nu2_hat*5, length = 100)
  loglikes = -sapply(log(nu2_vec), log_like, 
                     sigma2_vec = sigma2_vec, 
                     theta_hat_vec = theta_hat_vec)
  l_max = - log_like(sigma2_vec, log(nu2_hat), theta_hat_vec)
  likeli_df = data.frame(nu2 = nu2_vec, loglike = loglikes)
  
  ggplot(likeli_df, aes(nu2, loglike)) + 
    geom_line(size = 1, col = 'dark red') + 
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




Borrow_vs_No_Borrow_plot = function(output_object){
  Subgroup = paste0("Trt. ", 1:length(unique(output_object$CIs$Effect)))
  df_temp = data.frame(Center = rep(theta_vec, 2), Effect = rep(Subgroup, 2), 
                       Estimate = c(rep("Informed", K), rep("Uninformed", K)))
  cols = brewer.pal(n = 3, name = "Set1")
  colors = c("Informed" = cols[1], "Uninformed" = cols[2], "True" = cols[3])
  
  CI = output_object$CIs %>%
    mutate(Effect = plyr::mapvalues(Effect, 
                                    unique(output_object$CIs$Effect), 
                                    Subgroup),
           Estimate = as.character(Estimate))
  CI$Estimate[as.character(CI$Estimate) == 'Adjusted'] = 'Informed'
  CI$Estimate[as.character(CI$Estimate) == 'Unadjusted'] = 'Uninformed'
  
  meta_line = sum(theta_hat_vec/sigma2_vec)/sum(1/sigma2_vec)
  
  p = ggplot(CI, aes(Center, Effect, by = Estimate)) + 
    geom_errorbarh(mapping = aes(xmin = LL, xmax = UL, col = Estimate), 
                   height = .2, 
                   size = .5,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 3, mapping = aes(Center, Effect, col = Estimate),
               position = position_dodge(width = 0.5)) +
    xlab("Treatment effect") + 
    ylab("") + 
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          legend.position = "top",
          legend.text = element_text(size = 10), 
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(hjust = 1, size = 10, angle = 30),
          plot.title = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 12, face = "bold"),
          axis.line=element_line(),
          plot.margin = margin(0, 3, 0, 0, "cm")) + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_vline(xintercept = meta_line, 
               linetype = 'dashed', col = 'dark grey') + 
    geom_point(df_temp, mapping = aes(Center, Effect, col = "True"),
               position = position_dodge(width = 0.5), size = 2) + 
    scale_color_manual(name = "Trt. Effect", values = colors)
  
  return(p)
}




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
    scale_x_discrete(labels = TeX(unique(df$`Trt. Eff.`))) + 
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




MSE_vec = function(vec, sigma2_vec, theta_vec){
  MSE(vec[1], sigma2_vec, theta_vec, vec[2])
}




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




trt_eff_to_sd = function(n, p, theta){
  sqrt(1/n*(1/p*(1 + exp(-theta)) - 2))
}



trt_eff_to_sd_OR = function(n, p, theta){
  u = p*exp(theta)
  sig2 = 1/n*((2-p) + 1/u*(1 - p + u)^2)/(1 - p)
  
  sqrt(sig2)
}



trt_eff_to_sd_vec = function(n, vec){
  p = vec[1]
  theta = vec[2]
  sqrt(1/n*(1/p*(1 + exp(-theta)) - 2))
}



a_coeffs = function(theta_hat_vec, 
                    sigma2_vec, nu2, i){
  k = length(sigma2_vec)
  v = 1/(sigma2_vec + nu2)
  sum_v = sum(v)
  a = sigma2_vec[i]*v[i]*v
  a[i] = a[i] + nu2*v[i]*sum_v
  
  return(a/sum_v)
}




power_func = function(n, CER, RR, N_MC, theta_hat_vec, sigma2_vec, alpha){
  theta_new = as.numeric(log(RR))
  sigma2_new = as.numeric((trt_eff_to_sd(n, CER, theta_new))^2)
  
  CER_vec = qnorm(randtoolbox::sobol(N_MC), CER, sqrt(CER*(1 - CER)/n))
  theta_hat_new = qnorm(randtoolbox::sobol(N_MC), theta_new, sqrt(sigma2_new))
  sigma2_hat_new = (apply(cbind(CER_vec, theta_hat_new), 1,
                          trt_eff_to_sd_vec, n = n))^2
  
  df = cbind(theta_hat_new, sigma2_hat_new)
  
  power_aux_func = function(v){
    theta_hat_last = v[1]
    sigma2_hat_last = v[2]
    
    K = length(theta_hat_vec)
    
    sigma2_vec_aug = c(sigma2_vec, sigma2_hat_last)
    theta_vec = c(theta_hat_vec, theta_new)
    theta_hat_vec_aug = c(theta_hat_vec, theta_hat_last)
    
    nu2_hat = nu2_Emp_Bayes(sigma2_vec_aug, theta_hat_vec_aug)
    a = a_coeffs(theta_hat_vec_aug, sigma2_vec_aug, nu2_hat, K + 1)
    sigmastarKplus1 = sqrt(adjusted_variance(sigma2_vec_aug, nu2_hat, K+1))
    
    A = (sigmastarKplus1*qnorm(1 - alpha) + t(a)%*%theta_vec)/(a[K+1]*sqrt(sigma2_hat_last))
    
    return(c(1 - pnorm(A)))
  }
  
  
  pow = apply(df, 1, power_aux_func)
  return(mean(pow))
}





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






power_func_new = function(n, CER, RR, N_MC, theta_hat_vec, sigma2_vec, alpha){
  set.seed(n*N_MC)
  
  theta_new = as.numeric(log(RR))
  sigma2_new = as.numeric((trt_eff_to_sd(n, CER, theta_new))^2)
  
  CER_vec = qnorm(randtoolbox::sobol(N_MC), CER, sqrt(CER*(1 - CER)/n))
  theta_hat_new = qnorm(randtoolbox::sobol(N_MC), theta_new, sqrt(sigma2_new))
  sigma2_hat_new = (apply(cbind(CER_vec, theta_hat_new), 1,
                          trt_eff_to_sd_vec, n = n))^2
  
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






power_func_fast = function(theta_hat_vec, sigma2_vec,
                           n, CER, RR, alpha){
  K = length(theta_hat_vec)
  p_ctrl = CER
  theta = log(RR)
  
  sigmaKplus1 = trt_eff_to_sd(n, p_ctrl, theta)
  theta_hat_vec = c(theta_hat_vec, theta)
  sigma2_vec = c(sigma2_vec, sigmaKplus1^2)
  
  nu2_hat = nu2_Emp_Bayes(sigma2_vec, theta_hat_vec)
  a = a_coeffs(theta_hat_vec, sigma2_vec, nu2_hat, K + 1)
  
  sigmastarKplus1 = sqrt(adjusted_variance(sigma2_vec, nu2_hat, K+1))
  
  A = (sigmastarKplus1*qnorm(1 - alpha) + 
         t(a)%*%theta_hat_vec)/(a[K+1]*sigmaKplus1)
  
  return(c(1 - pnorm(A)))
}






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




beta_single_simulation = function(CER, RR, n, N_MC, alpha, sim_num){
  set.seed(sim_num)
  
  n_events = rbinom(2, n, CER*c(1, RR))
  post_eff = beta_post_efficacy(n_events, rep(n, 2), N_MC)
  
  return(post_eff > 1 - alpha)
}




beta_power_func = function(CER, RR, n, N_MC, alpha, N_sim){
  out = sapply(1:N_sim, 
               beta_single_simulation,
               CER = CER, RR = RR, n = n, 
               N_MC = N_MC, alpha = alpha)
  
  mean(out)
}




power_vs_n_df = function(n_vec, CER, RR, N_MC, theta_hat_vec, 
                         sigma2_vec, alpha, N_sim){
  cl = makeSOCKcluster(parallel::detectCores() - 1)
  registerDoSNOW(cl)
  
  pb = txtProgressBar(max = length(n_vec), style=3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress=progress)
  
  l = foreach(i = 1:length(n_vec),
              .options.snow = opts,
              .export = c(
                'nu2_Emp_Bayes',
                "trt_eff_to_sd",
                'trt_eff_to_sd_vec',
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
                'trt_eff_to_sd_vec',
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
                'trt_eff_to_sd_vec',
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






power_vs_RRR_plot = function(df, power_thresh){
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




weights_diff = function(sigma2_vec, nu2, target_weight){
  weights(sigma2_vec, nu2, 3)[1:2] - target_weight
}


weight_to_nu2 = function(sigma2_vec, target_weight){
  nu2_seq = seq(1e-6, .1, by = 1e-6)
  temp = t(sapply(nu2_seq, weights_diff, sigma2_vec = sigma2_vec, 
                  target_weight = target_weight))
  
  nu2_target = nu2_seq[max(apply(abs(temp), 2, which.min))]
  return(nu2_target)
}


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





library(ggplot2)
library(lemon)
library(coda)
library(mcmcplots)
library(dplyr)
library(rstan)
library(shinystan)
library(xtable)
library(scales)
library(RColorBrewer)
library(metafor)
library(cowplot)
library(kableExtra)
library(reshape2)
library(gtsummary)
library(LaplacesDemon)
library(doParallel)
library(doSNOW)
library(parallel)
library(tictoc)
library(OpenRepGrid)
library(randtoolbox)
library(grDevices)
library(latex2exp)
library(patchwork)
library(cobs)
library(tictoc)





theta_hat_vec = c(log(1.08) + log(2.19/1.08)/2, 
                  log(1.3) + log(2.59/1.3)/2, 
                  log(0.64) + log(3.46/0.64)/2)
sigma2_vec = c((log(2.19/1.08)/3.92)^2, (log(2.59/1.3)/3.92)^2, 
               (log(3.46/0.64)/3.92)^2)







weight_vec = seq(0, 1, by = .005)
out95 = t(sapply(weight_vec, output_fixed_w,
                 theta_hat_vec = theta_hat_vec,
                 sigma2_vec = sigma2_vec,
                 alpha = .05))

out98 = t(sapply(weight_vec, output_fixed_w,
                 theta_hat_vec = theta_hat_vec,
                 sigma2_vec = sigma2_vec,
                 alpha = .02))

df = as.data.frame(rbind(out95, out98)) %>% 
  mutate(Weight = rep(weight_vec, 2),
         Coverage = rep(c('95%', '98%'), each = length(weight_vec))) %>% 
  setNames(c('Mean', 'LL', 'UL', 'Weight', 'Coverage'))


min_weight95 = df$Weight[min(which(df$Coverage == '95%' & df$LL >= 0))]
min_weight98 = df$Weight[min(which(df$Coverage == '98%' & df$LL >= 0))]


post_prob_effic95 = df %>% 
  filter(Weight == .26 & Coverage == '95%') %>% 
  mutate(post_effic = pnorm(2*qnorm(.975)*Mean/(UL-LL)))

post_prob_effic95

post_prob_effic98 = df %>% 
  filter(Weight == Weight[which.min(abs(Weight - .35))] 
         & Coverage == '98%') %>% 
  mutate(post_effic = pnorm(2*qnorm(.99)*Mean/(UL-LL)))

post_prob_effic98


p = ggplot(df) + 
  geom_line(aes(Weight, LL, col = Coverage), 
            linetype = 'dashed', linewidth = 1) + 
  geom_line(aes(Weight, UL, col = Coverage), 
            linetype = 'dashed', linewidth = 1) + 
  geom_line(aes(Weight, Mean), 
            linetype = 'dashed', linewidth = 1) + 
  geom_hline(yintercept = 0, col = 'black') + 
  geom_segment(aes(x = min_weight95, xend = min_weight95,
                   y = -Inf, yend = 0)) + 
  geom_segment(aes(x = min_weight98, xend = min_weight98,
                   y = -Inf, yend = 0)) + 
  geom_point(aes(x = min_weight95, y = 0), 
             shape = 21, size = 2, fill = "white", 
             colour = "black", stroke = 1.5) + 
  geom_point(aes(x = min_weight98, y = 0), 
             shape = 21, size = 2, fill = "white", 
             colour = "black", stroke = 1.5) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        legend.position = 'top',
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=10),
        axis.line=element_line()) + 
  xlab(expression(omega)) + 
  ylab('log(OR)') + 
  scale_color_brewer(palette = 'Set1')

p



(nu2_95 = weight_to_nu2(sigma2_vec, min_weight95))
(w = weights(sigma2_vec, nu2_95, 3)[-3])
(ESS_increment95 = 93*sum(t(sigma2_vec[3]/sigma2_vec[-3])%*%w))


(nu2_98 = weight_to_nu2(sigma2_vec, min_weight98))
(w = weights(sigma2_vec, nu2_98, 3)[-3])
(ESS_increment98 = 93*sum(t(sigma2_vec[3]/sigma2_vec[-3])%*%w))





theta_hat_vec = c(log(1.08) + log(2.19/1.08)/2, 
                  log(1.3) + log(2.59/1.3)/2)
sigma2_vec = c((log(2.19/1.08)/3.92)^2, 
               (log(2.59/1.3)/3.92)^2)


n = 395
CER = .4
alpha = 1.03e-5
N_MC = 1e5

power_func_sim_OR(n, CER, OR = 1, N_MC, 
                  theta_hat_vec, sigma2_vec, alpha)

power_func_sim_OR(n, CER, OR = 1.5, N_MC, 
                  theta_hat_vec, sigma2_vec, alpha)


n = 50

weight_seq = seq(.15, .35, by = .005)
nu2_vec = sapply(weight_seq, weight_to_nu2,
                 sigma2_vec = c(sigma2_vec, 
                                (trt_eff_to_sd_OR(n, .4, log(1.5)))^2)
)

power975 = sapply(nu2_vec, power_func_sim_OR_fixed_nu2,
                  n = n, CER = CER, OR = 1.5, N_MC = 1e5,
                  theta_hat_vec = theta_hat_vec,
                  sigma2_vec = sigma2_vec, alpha = .025)

power99 = sapply(nu2_vec, power_func_sim_OR_fixed_nu2,
                 n = n, CER = CER, OR = 1.5, N_MC = 1e5,
                 theta_hat_vec = theta_hat_vec,
                 sigma2_vec = sigma2_vec, alpha = .01)


df = data.frame(weight = rep(weight_seq, 2),
                Power = c(power975, power99),
                Threshold = c(rep('97.5%', length(weight_seq)),
                              rep('99.0%', length(weight_seq))))

target_weight = df %>% 
  group_by(Threshold) %>% 
  summarize(target_weight = weight_seq[min(which(Power >= .8))]) %>% 
  as.data.frame() %>% 
  dplyr::select(target_weight) %>% 
  as.vector() %>% 
  unlist()

target_power = rbind(
  df %>% 
    filter(Threshold == "97.5%" & weight == target_weight[1]),
  df %>% 
    filter(Threshold == "99.0%" & weight == target_weight[2])
) %>% 
  dplyr::select(Power)

targets = data.frame(Threshold = c("97.5%", "99.0%"),
                     Power = target_power,
                     Weight = target_weight)

p = ggplot(df, aes(weight, Power, col = Threshold)) +
  geom_line(linewidth = 1) + 
  geom_segment(targets,
               mapping = aes(x = Weight, xend = Weight,
                             y = 0, yend = Power),
               linetype = 'dashed') +
  geom_segment(targets,
               mapping = aes(x = min(df$weight), xend = Weight,
                             y = Power, yend = Power),
               linetype = 'dashed') +
  geom_point(targets,
             mapping = aes(x = Weight, y = Power), 
             shape = 21, size = 2, fill = "white", stroke = 1.5) + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        legend.position = 'top',
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=10),
        axis.line=element_line()) + 
  xlab(expression(omega)) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_color_brewer(palette = "Set1")

p
