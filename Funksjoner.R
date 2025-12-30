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

