## Generate Data ========================================================================
gen_data <- function(treats, f_fun, sig, n_trials) {
  treats <- treats %>%
    mutate(m1 = f_fun(Ls-Lt1),
           m2 = f_fun(Ls-Lt2),
           dm = m1 - m2,
           
           pt1 = pnorm(dm, mean = 0, sd = sig))
  
  sim_data <- data.frame()
  for (i in 1:n_trials) {
    tmp <- treats
    tmp$rand <- runif(nrow(tmp), 0, 1)
    tmp$R <- 1
    tmp[which(tmp$rand < tmp$pt1),'R'] <- 0

    tmp <- tmp %>%
      dplyr::select(Ls, Lt1, Lt2, R)
    
    sim_data <- bind_rows(sim_data, tmp)
  }
  
  sim_data
  
}


## Format for Discrete MLE ==============================================================
format_data_discreteMLE <- function(data){
  data %>% 
    mutate(d1 = abs(Ls - Lt1),
           d2 = abs(Ls - Lt2),
           d1 = as.factor(d1),
           d2 = as.factor(d2))
}


## Discrete Likelihood with SD Free Parameter ===========================================
discrete.lklh.freeSD <- function(data, par){
  my_par <- c(0, cumsum(abs(par[1:(length(par)-1)])))
  my_par <- my_par / max(my_par)
  my_par <- my_par * FixedMax
    #if(max(my_par) > FixedMax){my_par <- FixedMax * my_par/max(my_par)}
  my_par <- c(my_par, FixedMax)
  my_sd <- par[length(par)]

  data %>%
    mutate(m1 = my_par[as.numeric(d1)],                       # grab the d1th parameter
           m2 = my_par[as.numeric(d2)],                       # grab the d2th parameter
           dm = m1 - m2,                                   # calculate diff
           
           pt1 = pnorm(dm, mean = 0, sd = my_sd),          # prob of test 1
           pt2 = 1 - pt1,                                  # prob of test 2
           
           pt1 = replace(pt1, pt1 == 0, .Machine$double.eps),      # correcting 0
           pt1 = replace(pt1, pt1 == 1, 1 - .Machine$double.eps),  # correcting 1
           pt2 = replace(pt2, pt2 == 0, .Machine$double.eps),      # correcting 0
           pt2 = replace(pt2, pt2 == 1, 1 - .Machine$double.eps),  # correcting 1
           
           logL = -log((1-R)*pt1 + R*pt2)) %>%             # calculating the NLL
    
    dplyr::select(logL) %>% sum()                          # summing the NLL
}

## Discrete Likelihood with Max Free Parameter ==========================================
discrete.lklh.freeMax <- function(data, par){
  #my_par <- c(0, cumsum(abs(par[1:(length(par)-1)])))
  my_par <- c(0, cumsum(abs(par)))
  my_sd <- FixedSDHat

  data %>%
    mutate(m1 = my_par[as.numeric(d1)],                       # grab the d1th parameter
           m2 = my_par[as.numeric(d2)],                       # grab the d2th parameter
           dm = m1 - m2,                                   # calculate diff
           
           pt1 = pnorm(dm, mean = 0, sd = my_sd),          # prob of test 1
           pt2 = 1 - pt1,                                  # prob of test 2
           
           pt1 = replace(pt1, pt1 == 0, .Machine$double.eps),      # correcting 0
           pt1 = replace(pt1, pt1 == 1, 1 - .Machine$double.eps),  # correcting 1
           pt2 = replace(pt2, pt2 == 0, .Machine$double.eps),      # correcting 0
           pt2 = replace(pt2, pt2 == 1, 1 - .Machine$double.eps),  # correcting 1
           
           logL = -log((1-R)*pt1 + R*pt2)) %>%             # calculating the NLL
    
    dplyr::select(logL) %>% sum()                          # summing the NLL
}
