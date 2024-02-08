#*******************************************************************************
#*******************************************************************************
#* The code used for the paper "Alone, together: on the benefits of Bayesian 
#* borrowing in a meta-analytic setting", to appear on Pharmaceutical Statistics 
#* (2023).
#* 
#* All code was written by Ofir Harari
#*******************************************************************************
#*******************************************************************************


library(ggplot2)
library(lemon)
library(coda)
library(mcmcplots)
library(dplyr)
library(rstan)
library(scales)
library(kableExtra)
library(reshape2)
library(LaplacesDemon)
library(doParallel)
library(doSNOW)
library(parallel)
library(latex2exp)
library(patchwork)
library(cobs)


wd = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
#Sourcing functions
source('Alone_Together_Functions.R')

################################################################################
################################################################################
## Part 1: posterior inference
################################################################################
################################################################################

#*******************************************************************************
#* Simulated data: used throughout the paper
#*******************************************************************************
seed = 20
set.seed(seed)

alpha = .05
K = 5
p_ctrl = .4
p_trt = .3
sigma2_ctrl = 0.2 
sigma2_trt = 0.2
n_ctrl = 75
n_trt = 75

data = random_data(K = K, p_ctrl = p_ctrl, p_trt = p_trt, 
                   sigma2_ctrl = sigma2_ctrl, 
                   sigma2_trt = sigma2_trt, 
                   n_ctrl = n_ctrl, n_trt = n_trt)

display_simulated_data(data)

theta_vec = data$theta_vec
theta_hat_vec = data$theta_hat_vec
sigma2_vec = data$sigma2_vec



#*******************************************************************************
#* Plotting the Bayesian borrowing fraction of studies 1-4 for the informing of
#* study 5 as a function of the study-to-study variance
#*******************************************************************************
trial_number = 5

# the study-to-study variance corresponding to a maximum borrowing fraction of 0.25
target_nu2 = .036 

# grid of variances for borrowing fraction evaluation
nu2_min = .01
nu2_max = .1
nu2_delta = .0001

weight_vs_nu2_plot(sigma2_vec, nu2_min, nu2_max, 
                   nu2_delta, trial_number, target_nu2)


#*******************************************************************************
#* Joint Bayesian analysis with a borrowing fraction not exceeding 0.25
#*******************************************************************************
out_max_BFP = output_fixed_nu(theta_hat_vec, sigma2_vec, 
                              alpha = .05, nu2 = target_nu2)


#*******************************************************************************
#* Joint Bayesian analysis using empirical Bayes to estimate  the study-to-study 
#* variance
#*******************************************************************************
out_Emp_Bayes = output(theta_hat_vec, sigma2_vec, alpha = .05)

nu2_hat = out_Emp_Bayes$nu2_hat
nu2_est_plot(sigma2_vec, theta_hat_vec, nu2_hat)


#*******************************************************************************
#* Full joint Bayesian analysis using a Stan program
#* If the code fails to run due to Stan version mismatch, visit - 
#* https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows#r-42
#* for instructions.
#*******************************************************************************
fit = stan_run(theta_hat_vec, sigma2_vec, N_iter = 1e5,  
               N_chains = 4, thin_fact = 10)

# Rhat and effective sample size diagnostics
summary(fit)$summary %>% 
  apply(., 2, round, digits = 2)

# Posterior nu2 histogram
nu2_posterior_plot(fit, x_lim = .3)

# Visual diagnostics
stan_diagnostics(fit)



#*******************************************************************************
#* Comparing the result of the different methods
#*******************************************************************************
Comparison = Method_comparison(fit, out_Emp_Bayes, theta_vec, 
                               theta_hat_vec, out_max_BFP)

Comparison$Table %>% 
  kable(align = c('c', rep('r', ncol(Comparison$Table) - 1)), 
        'html', escape = F) %>%
  kable_styling(full_width = F, font_size = 14,
                bootstrap_options = c("striped")) %>%
  row_spec(0, font_size = 14, align = 'c')


Comparison$Plot



################################################################################
################################################################################
## Part 2: estimation analysis
################################################################################
################################################################################

#*******************************************************************************
#* Calculating and plotting theoretical RMSEs for all parameters 
#*******************************************************************************
nu2_vec = seq(0, .75, length = 5000)
p = RMSE_vs_nu2_plot_table(nu2_vec, sigma2_vec, nu2_hat, 
                           theta_vec, theta_hat_vec)
p$Plot
p$Table


#*******************************************************************************
#* The simulation from Example 3.2: simulating with low and high heterogeneity
#* and summarizing the results
#*******************************************************************************
sim_small_nu2 = adj_vs_unadj_simulation(K = 5, 
                                        p_ctrl = .4, 
                                        p_trt = .3, 
                                        sigma2_ctrl = .2, 
                                        sigma2_trt = .2, 
                                        n_ctrl = 75, 
                                        n_trt = 75, 
                                        alpha = .05, 
                                        N_sims = 1e5) %>% 
  mutate(Heterogeneity = 'Low')


sim_large_nu2 = adj_vs_unadj_simulation(K = 5, 
                                        p_ctrl = .4, 
                                        p_trt = .3, 
                                        sigma2_ctrl = .4, 
                                        sigma2_trt = .4, 
                                        n_ctrl = 75, 
                                        n_trt = 75, 
                                        alpha = .05, 
                                        N_sims = 1e5) %>% 
  mutate(Heterogeneity = 'High')

sim = rbind(sim_small_nu2, sim_large_nu2)
sim[sim == Inf] = NA

sim %>% 
  mutate(Heterogeneity = factor(Heterogeneity, 
                                levels = c('Low', 'High'))) %>% 
  group_by(Heterogeneity, Type) %>% 
  summarize(Bias = mean(delta, na.rm = T),
            RMSE = sqrt(mean(delta^2, na.rm = T)),
            Length = mean(CI_length, na.rm = T),
            Coverage = mean(covered, na.rm = T)) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  mutate(Coverage = paste0(100*Coverage, '%')) %>% 
  kable(align = 'c') %>%
  kable_styling(full_width = F, font_size = 14,
                bootstrap_options = c("striped")) %>%
  row_spec(0, font_size = 14) 




################################################################################
## Design (power and sample size calculations)
################################################################################


#*******************************************************************************
#* Type I error rate and power for the informed model
#*******************************************************************************
power_func_sim(n = 339, CER = .4, RR = 1, N_MC = 1e5, 
               theta_hat_vec, sigma2_vec, alpha = .008090097)
# 2.5% type I error rate for the informed model

power_func_sim(n = 339, CER = .4, RR = .75, N_MC = 1e5, 
               theta_hat_vec, sigma2_vec, alpha = .008090097)
# 80% power for the informed model


#*******************************************************************************
#* Type I error rate and power for the uninformed model
#*******************************************************************************
beta_power_func(CER = .4, RR = 1, n = 478, N_MC = 5e3, 
                alpha = .008090097, N_sim = 1e5)
# 0.865% type I error rate for the uninformed model
beta_power_func(CER = .4, RR = .75, n = 478, N_MC = 5e3, 
                alpha = .008090097, N_sim = 1e5)
# 80% power for the uninformed model


#*******************************************************************************
#* Type I error rate and power for the uninformed model after re-calibrating to 
#* 2.5% type I error rate
#*******************************************************************************
beta_power_func(CER = .4, RR = 1, n = 360, N_MC = 5e3, 
                alpha = .0245, N_sim = 5e4)
# 2.5% type I error rate for the uninformed model
beta_power_func(CER = .4, RR = .75, n = 360, N_MC = 5e3, 
                alpha = .0245, N_sim = 5e4)
# 80% power for the uninformed model




#*******************************************************************************
#* Comparing power (vs. sample size) between informed and uninformed (and
#* non-calibrated) Bayesian models
#*******************************************************************************
RR = .75
CER = .4
N_MC = 2.5e3
alpha = .008090097
N_sim = 3e4
n_vec = seq(199, 750, by = 4)

monkey = 0 
# do not change "monkey" to 1 unless you want to re-run the (extremely long)
# power simulation

if(monkey){
  df = power_vs_n_df(n_vec, CER, RR, N_MC, theta_hat_vec, 
                     sigma2_vec, alpha, N_sim)
  
  save(df, file = "power_samp_size_df.Rda")
}


load("power_samp_size_df.Rda") # reading a pre-calculated data frame
power_vs_n_plot(df, power_thresh = 0.8)

# The model with borrowing will require (in this example) 339 patients per arm
# for 80% power, compared to 478 per arm by the standard beta-binomial model
# without borrowing.


#*******************************************************************************
#* Comparing power (vs. relative risk reduction) between informed and uninformed 
#* (and non-calibrated) Bayesian models
#*******************************************************************************
RR_vec = seq(.675, .95, length = 200) # we let the relative risk vary
n = 500 # we fix the per-arm sample size on 500 patients

monkey = 0
# do not change "monkey" to 1 unless you want to re-run the (extremely long)
# power simulation

if(monkey){
  df = power_vs_RRR_df(n, CER, RR_vec, N_MC, theta_hat_vec, 
                       sigma2_vec, alpha, N_sim)
  
  save(df, file = "power_rel_risk_reduct_df.Rda")
}

load("power_rel_risk_reduct_df.Rda") # reading a pre-calculated data frame
power_vs_RRR_plot(df, power_thresh = 0.8)

# The model with borrowing will have 80% power to detect an RRR of 21.7%,
# compared to 24.5% by the standard beta-binomial model without borrowing.



################################################################################
## Systemic Lupus Erythematosus example: design
################################################################################

#*******************************************************************************
#* Log odds-ratios of NCT00410384, NCT00424476 and NCT01649765, by order -
#*******************************************************************************
theta_hat_vec = c(log(1.08) + log(2.19/1.08)/2, 
                  log(1.3) + log(2.59/1.3)/2, 
                  log(0.64) + log(3.46/0.64)/2)

#*******************************************************************************
#* The corresponding variances, derived from the published confidence intervals
#*******************************************************************************
sigma2_vec = c((log(2.19/1.08)/3.92)^2, 
               (log(2.59/1.3)/3.92)^2, 
               (log(3.46/0.64)/3.92)^2)

#*******************************************************************************
#* Calibrating to type I error rate constraints
#*******************************************************************************  
power_func_sim_OR(n = 396, CER = .4, OR = 1, N_MC = 1e5, 
                  theta_hat_vec, sigma2_vec, alpha = 1.05e-5)

power_func_sim_OR(n = 396, CER = .4, OR = 1.5, N_MC = 1e5, 
                  theta_hat_vec, sigma2_vec, alpha = 1.05e-5)


#*******************************************************************************
#* Fixing the sample size and calculating the between-study variances that 
#* correspond to a vector of borrowing fractions
#******************************************************************************* 
n = 50 # number of patients per arm

weight_seq = seq(.15, .35, by = .005)

nu2_vec = sapply(weight_seq, weight_to_nu2,
                 sigma2_vec = c(sigma2_vec, 
                                (trt_eff_to_sd_OR(n, .4, log(1.5)))^2)
)



monkey = 0 
# do not change "monkey" to 1 unless you want to re-run the (extremely long)
# power simulation

if(monkey){
#*******************************************************************************
#* calculating the power with respect to 97.5% and 99% posterior probability of 
#* efficacy for every such borrowing fraction
#******************************************************************************* 
  power975 = sapply(nu2_vec, power_func_sim_OR_fixed_nu2,
                    n = n, CER = .4, OR = 1.5, N_MC = 1e5,
                    theta_hat_vec = theta_hat_vec,
                    sigma2_vec = sigma2_vec, alpha = .025)
  
  power99 = sapply(nu2_vec, power_func_sim_OR_fixed_nu2,
                   n = n, CER = .4, OR = 1.5, N_MC = 1e5,
                   theta_hat_vec = theta_hat_vec,
                   sigma2_vec = sigma2_vec, alpha = .01)
  
  
#*******************************************************************************
#* storing the results in a data frame and find the minimum borrowing fractions 
#* that lead to meeting the thresholds
#******************************************************************************* 
  df = data.frame(weight = rep(weight_seq, 2),
                  Power = c(power975, power99),
                  Threshold = c(rep('97.5%', length(weight_seq)),
                                rep('99.0%', length(weight_seq))))
  
  save(df, file = "Lupus_power_df.Rda")
}

load("Lupus_power_df.Rda") # reading a pre-calculated data frame

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

(targets = data.frame(Threshold = c("97.5%", "99.0%"),
                      Power = target_power,
                      Weight = target_weight))


#*******************************************************************************
#* Plotting the power curves
#******************************************************************************* 
ggplot(df, aes(weight, Power, col = Threshold)) +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::percent) +
  scale_color_brewer(palette = "Set1")



################################################################################
################################################################################
## Part 3: Systemic Lupus Erythematosus example: analysis
################################################################################
################################################################################

#*******************************************************************************
#* jointly analyzing the three trials for a vector of borrowing fractions, 
#* calculating both 95% and 98% credible intervals
#******************************************************************************* 
weight_vec = seq(0, 1, by = .005)

out95 = t(sapply(weight_vec, output_fixed_w,
                 theta_hat_vec = theta_hat_vec,
                 sigma2_vec = sigma2_vec,
                 alpha = .05))

out98 = t(sapply(weight_vec, output_fixed_w,
                 theta_hat_vec = theta_hat_vec,
                 sigma2_vec = sigma2_vec,
                 alpha = .02))

#*******************************************************************************
#* storing the results in a data frame and calculating the minimum borrowing 
#* fractions required for the intervals to be entirely to the right of the null, 
#* as well as the posterior probabilities of efficacy
#******************************************************************************* 
df = as.data.frame(rbind(out95, out98)) %>% 
  mutate(Weight = rep(weight_vec, 2),
         Coverage = rep(c('95%', '98%'), each = length(weight_vec))) %>% 
  setNames(c('Mean', 'LL', 'UL', 'Weight', 'Coverage'))

# Is this smaller than 0.26?
(min_weight95 = df$Weight[min(which(df$Coverage == '95%' & df$LL >= 0))])
#Is this smaller than 0.35?
(min_weight98 = df$Weight[min(which(df$Coverage == '98%' & df$LL >= 0))])


#*******************************************************************************
#* Plotting the results
#*******************************************************************************
ggplot(df) + 
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
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)),
                     labels = percent) + 
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


#*******************************************************************************
#* Calculating the posterior probability of efficacy at the borrowing fraction 
#* we preset at the design stage, and confirming that it crosses the 97.5% (99%) 
#* threshold
#*******************************************************************************
post_prob_effic95 = df %>% 
  filter(Weight == .26 & Coverage == '95%') %>% 
  mutate(post_effic = pnorm(2*qnorm(.975)*Mean/(UL-LL)))

post_prob_effic95

post_prob_effic98 = df %>% 
  filter(Weight == Weight[which.min(abs(Weight - .35))] 
         & Coverage == '98%') %>% 
  mutate(post_effic = pnorm(2*qnorm(.99)*Mean/(UL-LL)))

post_prob_effic98


#*******************************************************************************
#* Calculating the effective number of patients borrowed from the adult trials
#* for a 97.5% posterior probability of efficacy threshold
#*******************************************************************************
nu2_95 = weight_to_nu2(sigma2_vec, min_weight95)
(w = weights(sigma2_vec, nu2_95, 3)[-3])
(ESS_increment95 = 93*sum(t(sigma2_vec[3]/sigma2_vec[-3])%*%w))

#*******************************************************************************
#* Calculating the effective number of patients borrowed from the adult trials
#* for a 99% posterior probability of efficacy threshold
#*******************************************************************************
nu2_98 = weight_to_nu2(sigma2_vec, min_weight98)
(w = weights(sigma2_vec, nu2_98, 3)[-3])
(ESS_increment98 = 93*sum(t(sigma2_vec[3]/sigma2_vec[-3])%*%w))
