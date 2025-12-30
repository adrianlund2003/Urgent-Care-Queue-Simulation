CI = function(est,alpha){
  n = length(est)
  mu = mean(est)
  s = sd(est)
  q = qt(1 - (1-alpha)/2,n-1)
  return(c(mu-q*s/sqrt(n),mu+q*s/sqrt(n)))
}

UCC = function(lambda,mu){
  t = 0
  n_in_queue = 0
  queue_list = c(n_in_queue)
  t_list = c(t)
  t_diff = c(t)
  
  while(t <= (50 * 24)){
    
    if(n_in_queue == 0){
      arrival_time = rexp(1,lambda)
      t = t + arrival_time
      n_in_queue = n_in_queue + 1
      queue_list = c(queue_list,n_in_queue)
      t_list = c(t_list,t)
      t_diff = c(t_diff,arrival_time)
    }
    
    if(n_in_queue != 0){
      treatment_time = rexp(1,mu)
      t = t + treatment_time
      n_in_queue = n_in_queue + rpois(1,treatment_time * lambda) - 1
      queue_list = c(queue_list,n_in_queue)
      t_list = c(t_list,t)
      t_diff = c(t_diff,treatment_time)
    }
  }
  
  weighted_sum = 0
  for (i in 1:length(queue_list)){
    weighted_sum = weighted_sum + queue_list[i] * t_diff[i]
  }
  expected_number_in_queue = weighted_sum / t
  expected_time_in_queue = expected_number_in_queue / lambda
  
  return(list(expected_time_in_queue = expected_time_in_queue,t_list=t_list,queue_list = queue_list))
}


gaussian_process <- function(theta_grid, theta_obs, y_obs, nugget = 1e-10) {
  m0 <- 0.5
  sigma2 <- 0.5^2
  
  corr <- function(theta1, theta2) {
    d <- abs(outer(theta1, theta2, "-"))
    return ((1+15*d)*exp(-15*d))
  }
  
  mean_vec_grid <- rep(m0, length(theta_grid))
  mean_vec_obs <- rep(m0, length(theta_obs))
  
  cov_gg <- sigma2 * corr(theta_grid, theta_grid)
  cov_go <- sigma2 * corr(theta_grid, theta_obs)
  cov_oo <- sigma2 * corr(theta_obs, theta_obs)
  
  cov_oo <- cov_oo + nugget * diag(length(theta_obs))
  R <- chol(cov_oo)                         
  cov_oo_inv <- chol2inv(R)  
  
  post_mean <- mean_vec_grid + cov_go %*% cov_oo_inv %*% (y_obs - mean_vec_obs)
  post_cov <- cov_gg - cov_go %*% cov_oo_inv %*% t(cov_go)
  post_sd <- sqrt(diag(post_cov))
  
  alpha <- 0.1
  z <- qnorm(1-alpha/2)
  
  lower_bound <- post_mean - z * post_sd
  upper_bound <- post_mean + z * post_sd
  
  
  plot(theta_grid, post_mean, type = "n", lwd = 2,
       ylab = expression(Y(theta)), xlab = expression(theta),
       ylim = c(0, 1))
  
  polygon(c(theta_grid, rev(theta_grid)),
          c(lower_bound, rev(upper_bound)),
          col = rgb(0.2, 0.4, 1, 0.2), border = NA)
  
  lines(theta_grid, post_mean, lwd = 2)
  lines(theta_grid, lower_bound, lty = 2, col = "gray40")
  lines(theta_grid, upper_bound, lty = 2, col = "gray40")
  points(theta_obs, y_obs, pch = 19, col = "red")
  
  return(list(
    post_mean = post_mean,
    post_sd   = post_sd,
    post_cov  = post_cov
  ))
}

plot_prob <- function(threshold, post_mean, post_sd) {
  prob_below <- pnorm((threshold - post_mean) / post_sd)
  
  plot(theta_grid, prob_below, type = "l", lwd = 2, col = "blue",
       xlab = expression(theta),
       ylab = expression(P(Y(theta) < 0.30)),
       ylim = c(0, 1))
  abline(h = 0.5, lty = 2, col = "gray70")
  
  return (prob_below)
}


UCC_priority = function(lambda,mu,p){
  
  t = 0
  n_in_queue = 0
  u_in_queue = 0
  leftover_time = 0
  half_treated = FALSE
  
  n_in_queue_list = c(n_in_queue)
  u_in_queue_list = c(u_in_queue)
  t_list = c(t)
  t_diff = c(t)
  
  while (t <= 24*50){
    
    if (n_in_queue == 0 && u_in_queue == 0){
      arrival_time = rexp(1,lambda)
      dice = runif(1,0,1)
      if (dice <= p){
        u_in_queue = u_in_queue + 1
      } else {
        n_in_queue = n_in_queue + 1
      }
      t = t + arrival_time
      t_list = c(t_list,t)
      t_diff = c(t_diff,arrival_time)
      u_in_queue_list = c(u_in_queue_list,u_in_queue)
      n_in_queue_list = c(n_in_queue_list,n_in_queue)
    }
    
    if (u_in_queue != 0){
      
      treatment_time = rexp(1,mu)
      
      while(treatment_time > 0){
        arrival_time = rexp(1,lambda)
        if(treatment_time > arrival_time){
          dice = runif(1,0,1)
          if (dice <= p){
            u_in_queue = u_in_queue + 1
            t = t + arrival_time
            t_list = c(t_list,t)
            treatment_time = treatment_time - arrival_time
            t_diff = c(t_diff,arrival_time)
            u_in_queue_list = c(u_in_queue_list,u_in_queue)
            n_in_queue_list = c(n_in_queue_list,n_in_queue)
          } else {
            n_in_queue = n_in_queue + 1
            t = t + arrival_time
            treatment_time = treatment_time - arrival_time
            t_list = c(t_list,t)
            t_diff = c(t_diff,arrival_time)
            u_in_queue_list = c(u_in_queue_list,u_in_queue)
            n_in_queue_list = c(n_in_queue_list,n_in_queue)
          }
        } else {
          break
        }
      }
      
      u_in_queue = u_in_queue - 1
      t = t + treatment_time
      t_list = c(t_list,t)
      t_diff = c(t_diff,treatment_time)
      u_in_queue_list = c(u_in_queue_list,u_in_queue)
      n_in_queue_list = c(n_in_queue_list,n_in_queue)
    }
    
    if (u_in_queue == 0 && n_in_queue != 0){
      
      if (half_treated == TRUE){
        treatment_time = leftover_time
      } else {
        treatment_time = rexp(1,mu)
      }
      
      u_arrived = FALSE
      
      while(treatment_time > 0 && u_arrived == FALSE){
        arrival_time = rexp(1,lambda)
        if(treatment_time > arrival_time){
          dice = runif(1,0,1)
          if (dice <= p){
            u_in_queue = u_in_queue + 1
            treatment_time = treatment_time - arrival_time
            u_arrived = TRUE
          } else {
            n_in_queue = n_in_queue + 1
            treatment_time = treatment_time - arrival_time
            t = t + arrival_time
            t_list = c(t_list,t)
            t_diff = c(t_diff,arrival_time)
            u_in_queue_list = c(u_in_queue_list,u_in_queue)
            n_in_queue_list = c(n_in_queue_list,n_in_queue)
          }
        } else {
          break
        }
      }
      
      if (u_arrived == TRUE){
        leftover_time = treatment_time
        half_treated = TRUE
        t = t + treatment_time
        t_diff = c(t_diff,treatment_time)
        t_list = c(t_list,t)
        u_in_queue_list = c(u_in_queue_list,u_in_queue)
        n_in_queue_list = c(n_in_queue_list,n_in_queue)
        
      } else {
        n_in_queue = n_in_queue - 1
        half_treated = FALSE
        leftover_time = 0
        t = t + treatment_time
        t_list = c(t_list,t)
        t_diff = c(t_diff,treatment_time)
        u_in_queue_list = c(u_in_queue_list,u_in_queue)
        n_in_queue_list = c(n_in_queue_list,n_in_queue)
      }
    }
  }
  
  weighted_sum_u = 0
  weighted_sum_n = 0
  for (i in 1:length(u_in_queue_list)){
    weighted_sum_u = weighted_sum_u + u_in_queue_list[i] * t_diff[i]
    weighted_sum_n = weighted_sum_n + n_in_queue_list[i] * t_diff[i]
  }
  expected_u_in_queue = weighted_sum_u / t
  expected_n_in_queue = weighted_sum_n / t
  expected_w_u = expected_u_in_queue / (p * lambda)
  expected_w_n = expected_n_in_queue / ((1-p) * lambda)
  
  return(list(expected_w_u = expected_w_u,expected_w_n = expected_w_n,t_list = t_list,u_in_queue_list = u_in_queue_list,n_in_queue_list = n_in_queue_list))
}

