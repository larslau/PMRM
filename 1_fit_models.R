library(nlme)
library(contrast)

#
# ADNI decline and covariance (case study 1)
#

# Mean ADAS-cog scores 0, 6, 12, 18, 24, and 36 months after baseline
mean_pbo <- c(19.6, 20.5, 20.9, 22.7, 23.8, 27.4)

# Active arm (20% slowing of time-progression)
mean_act <- approx(x = c(0, 6, 12, 18, 24, 36),
                   y = mean_pbo,
                   xout = 0.8 * c(0, 6, 12, 18, 24, 36))$y

# Covariance matrix
cov_adni <- structure(c(45.15, 39.99, 45.1, 54.95, 53.58, 60.82, 
                        39.99, 57.78, 54.38, 66.33, 64.1, 74.67, 
                        45.1, 54.38, 72.01, 79.97, 77.64, 93.11, 
                        54.95, 66.33, 79.97, 109.77, 99.29, 121.66, 
                        53.58, 64.1, 77.64, 99.29, 111.41, 127.83, 
                        60.82, 74.67, 93.11, 121.66, 127.83, 191.41), .Dim = c(6L, 6L))
# Square root of covariance
cov_sqrt <- t(chol(cov_adni))


#
# Simulate trial with 500 subjects per arm
#

source('0_simulate_data.R')
set.seed(123)
dat <- simulate_trial(n_arm = 500,
                      mean_pbo = mean_pbo,
                      mean_act = mean_act,
                      cov_sqrt = cov_sqrt)

########################
# Fit different models #
########################

#
# MMRM/cLDA model
#

mmrm <- gls(model = y ~ 0 + act_vis,
            data = dat,
            correlation = corSymm(form = ~ as.numeric(visit) | id),
            weights = varIdent(form = ~ 1 | visit),
            method = 'ML')

# Test treatment differences
contrast(mmrm,
         a = list('act_vis' = 'act.36'),
         b = list('act_vis' = 'pbo.36'))

# #
# # Alternative MMRM/cLDA model with lme4 (faster for larger datasets)
# #
# 
# library(lme4)
# 
# mmrm_alt <- lmer(y ~ act_vis + 0 + (as.factor(visit) + 0 | id),
#              data = dat,
#              control = lmerControl(check.nobs.vs.nRE = 'ignore',
#                                    optimizer = 'optimx',
#                                    optCtrl = list(method = 'L-BFGS-B')))
# 
# # NOTE: lmer is overparametrized with a residual error term which
# # will likely give a warning about convergence that can safely be ignored
# #
# # However, some things are important to consider:
# #
# # * The model believes it has +1 variance parameter
# #   than what it effectively does and calculations that 
# #   rely on this (e.g. AIC, BIC, certain test statistics)
# #   should be computed or adjusted manually
# # * The residual variance parameter and the random effect
# #   variance-covariance matrix should never be interpreted 
# #   on their own, but once their combined effects have been
# #   added together
# 
# # Compare
# rbind(coef(mmrm), 
#       fixef(mmrm_alt))


#
# Proportional decline PMRM
#

# Mean function
PD_PMRM <- function(t, v0, v1, v2, v3, v4, v5, b) {
  months <- c(0, 6, 12, 18, 24, 36)
  v <- c(v0[1], v1[1], v2[1], v3[1], v4[1], v5[1])
  
  (1 - b) * (v[match(t, months)] - v[1]) + v[1]
}

# Initialize reduced decline to average observed
b_init <- 1 - mean(coef(mmrm)[c(2, 4, 6, 8, 10)] / coef(mmrm)[c(3, 5, 7, 9, 11)])

# Pull out placebo arm parameters
start_vec <- coef(mmrm)
start_vec <- c(start_vec[grepl('pbo', names(start_vec))], 
               b = b_init)

# Fit model
pd_pmrm <- gnls(model = y ~ PD_PMRM(M, v0, v1, v2, v3, v4, v5, b),
                data = dat,
                params = list(v0 + v1 + v2 + v3 + v4 + v5 ~ 1,
                              b ~ act + 0),
                correlation = corSymm(form = ~ as.numeric(visit) | id),
                weights = varIdent(form = ~ 1 | visit),
                start = start_vec,
                control = gnlsControl(nlsTol = 1))

# Test decline
anova(pd_pmrm)

#
# Time-PMRM
#

# Mean function
TPMRM <- function(t, v0, v1, v2, v3, v4, v5,
                  b1, b2, b3, b4, b5) {
  months <- c(0, 6, 12, 18, 24, 36)
  b <- cbind(0, b1, b2, b3, b4, b5) # b represents % slowing at each visit
  t_out <- (1 - b[cbind(1:length(t), match(t, months))]) * t 
  
  spline(x = months,
         y = c(v0[1], v1[1], v2[1], v3[1], v4[1], v5[1]),
         method = 'natural',
         xout = t_out)$y
}

# Pull out placebo arm parameters
start_vec <- coef(mmrm)
start_vec <- c(start_vec[grepl('pbo', names(start_vec))], 
               b = rep(0, 5))
# Fit model
t_pmrm <- gnls(model = y ~ TPMRM(M, v0, v1, v2, v3, v4, v5,
                                 b1, b2, b3, b4, b5),
               data = dat,
               params = list(v0 + v1 + v2 + v3 + v4 + v5 ~ 1,
                             b1 + b2 + b3 + b4 + b5 ~ act + 0),
               correlation = corSymm(form = ~ as.numeric(visit) | id),
               weights = varIdent(form = ~ 1 | visit),
               start = start_vec,
               control = gnlsControl(nlsTol = 1))

# Test slowing
anova(t_pmrm)


#
# Proportional slowing Time-PMRM
#

# Pull out placebo arm parameters
start_vec <- coef(mmrm)
start_vec <- c(start_vec[grepl('pbo', names(start_vec))], 
               b = 0)
# Fit model
pst_pmrm <- gnls(model = y ~ TPMRM(M, v0, v1, v2, v3, v4, v5,
                                   b, b, b, b, b),
                 data = dat,
                 params = list(v0 + v1 + v2 + v3 + v4 + v5 ~ 1,
                               b ~ act + 0),
                 correlation = corSymm(form = ~ as.numeric(visit) | id),
                 weights = varIdent(form = ~ 1 | visit),
                 start = start_vec,
                 control = gnlsControl(nlsTol = 1))

# Test slowing
anova(pst_pmrm)

# Testing proportional slowing assumption
anova(pst_pmrm, t_pmrm)

#
# Compare model fits
#

AIC(mmrm, 
    pd_pmrm, 
    t_pmrm, 
    pst_pmrm)



