#
# Code for simulating trial with visits at 0, 6, 12, 18, 24 months
#

simulate_trial <- function(n_arm = 500, 
                           M = c(0, 6, 12, 18, 24, 36), 
                           drop_out = 0,
                           mean_pbo, 
                           mean_act, 
                           cov_sqrt) {
  
  # Simulate data
  m <- length(M)
  y_pbo <- cov_sqrt %*% matrix(rnorm(m * n_arm), nrow = m) + mean_pbo
  y_act <- cov_sqrt %*% matrix(rnorm(m * n_arm), nrow = m) + mean_act
  
  dat <- data.frame(id = rep(1:(2 * n_arm), each = m), # Patient IDs
                    visit = factor(rep(1:m, 2 * n_arm)), # Visits
                    M = rep(M, n_arm), # Months since baseline
                    y = c(as.numeric(y_pbo), as.numeric(y_act)), # Outcome measure
                    trt = factor(rep(c('pbo', 'act'), each = n_arm * m)), # Treatment arm
                    act = rep(c(0, 1), each = n_arm * m)) # Dummy treatment arm
  dat$mod_trt <- dat$trt
  dat$mod_trt[dat$M == 0] <- 'pbo' # Cast 
  dat$act <- as.numeric(dat$trt == 'act')
  dat$act_vis <- with(dat, interaction(mod_trt, M))
  
  # Implement per-visit drop out
  for (v in M) {
    ids <- unique(subset(dat, M >= v)$id)
    # Sample 
    id_rm <- sample(ids, 
                    size = floor(length(ids) * drop_out))
    dat <- subset(dat, !(M >= v & (id %in% id_rm)))
  }
  
  return(dat) 
}
