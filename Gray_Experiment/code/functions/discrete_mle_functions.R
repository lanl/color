format_data_discreteMLE <- function(data){
  data %>% 
    mutate(d1 = abs(Ls - Lt1),
           d2 = abs(Ls - Lt2),
           d1 = as.factor(d1),
           d2 = as.factor(d2))
}

discrete.lklh <- function(data, par){
  my_par <- c(0, cumsum(abs(par[1:(length(par)-1)])))
    if(max(my_par) > 1){my_par <- my_par/max(my_par)}
  my_par <- c(my_par, 1)
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

multi_approx_dis <- function(data, n_approx, n_level){
  found_par_df <- data.frame()
  scale_par_df <- data.frame()
  sd_vec <- c()
  L_vec <- c()

  for(i in 1:n_approx) {
    init_pars <- c(runif(n_level-2, 0, 1/n_level), runif(1, 0.1, 1))
    test <- optim(par = init_pars,
                  fn = discrete.lklh,
                  data = data)

    L_vec <- c(L_vec, test$value)
    sd_vec <- c(sd_vec, test$par[length(test$par)])
        
    my_par <- c(0, cumsum(abs(test$par[1:(length(test$par)-1)])))
    if(max(my_par) > 1){my_par <- my_par/max(my_par)}
    my_par <- c(my_par, 1)
   
    found_row <- my_par
    scale_row <- rescale(found_row, to = c(0,1))
    print(found_row)
    
    found_par_df <- rbind(found_par_df, found_row)
    scale_par_df <- rbind(scale_par_df, scale_row)
  }
  my_list <- list('sd_vec' = sd_vec, 'L_vec' = L_vec, 
                  'found_par_df' = found_par_df,
                  'scale_par_df' = scale_par_df)

  my_list
}

visualizing_scaled_results <- function(data, scale_par_df, L_vec, n_levels){
  L_vec_min <- min(L_vec)
  L_vec_max <- max(L_vec)
  L_vec_ran <- L_vec_max - L_vec_min
  
  L_vec_cor <- -L_vec + L_vec_ran + 2*L_vec_min
  
  estimated_scaled <- sweep(scale_par_df, MARGIN = 1, L_vec_cor, '*')
  estimated <- colSums(estimated_scaled) / sum(L_vec)

  #estimated <- colMeans(estimated_scaled)
  CIlowerbd <- c()
  CIupperbd <- c()

  
  for(i in 2:(n_levels-1)){
    l <- t.test(scale_par_df[,i])$'conf.int'[1]
    u <- t.test(scale_par_df[,i])$'conf.int'[2]

    CIlowerbd <- c(CIlowerbd, l)
    CIupperbd <- c(CIupperbd, u)
  }

  CIlowerbd <- c(0, CIlowerbd, 1)
  CIupperbd <- c(0, CIupperbd, 1)
  print(CIlowerbd)
  print(CIupperbd)

  truelabel <- as.numeric(levels(data$d1))

  vis_df <- cbind(estimated, CIlowerbd, CIupperbd, truelabel)
  vis_df <- as.data.frame(vis_df)
  names(vis_df) <- c('estimated', 'CIlowerbd', 'CIupperbd', 'truelabel')
  
  raw_df <- as.data.frame(t(scale_par_df))
  names(raw_df) <- L_vec
  raw_df <- cbind(raw_df, truelabel)
  raw_df_m <- melt(raw_df, id.vars = 'truelabel')
  raw_df_m$variable <- as.numeric(as.character(raw_df_m$variable))
  
  p1 <- ggplot(vis_df, aes(x = truelabel, y = estimated)) + 
		geom_point(size = 4) + geom_errorbar(aes(ymin = CIlowerbd, ymax = CIupperbd), size = 1) +
                geom_line(data = raw_df_m, aes(x = truelabel, y = value, 
                                               group = as.factor(variable),
                                               color = variable)) +
                xlab('Difference in Absolute Units') + ylab('Estimated Perceived Difference') +
		scale_color_continuous(name = 'Negative Log\nLikelihood') #+
		#ggtitle('Discrete Estimation of Differences')

  p1

}

visualizing_agg_results <- function(data, scale_par_df, L_vec, n_levels){
  L_vec_min <- min(L_vec)
  L_vec_max <- max(L_vec)
  L_vec_ran <- L_vec_max - L_vec_min
  
  L_vec_cor <- -L_vec + L_vec_ran + 2*L_vec_min
  
  estimated_scaled <- sweep(scale_par_df, MARGIN = 1, L_vec_cor, '*')
  estimated <- colSums(estimated_scaled) / sum(L_vec)

  #estimated <- colMeans(estimated_scaled)
  CIlowerbd <- c()
  CIupperbd <- c()

  
  for(i in 2:(n_levels-1)){
    l <- t.test(scale_par_df[,i])$'conf.int'[1]
    u <- t.test(scale_par_df[,i])$'conf.int'[2]

    CIlowerbd <- c(CIlowerbd, l)
    CIupperbd <- c(CIupperbd, u)
  }

  CIlowerbd <- c(0, CIlowerbd, 1)
  CIupperbd <- c(0, CIupperbd, 1)

  truelabel <- as.numeric(levels(data$d1))

  vis_df <- cbind(estimated, CIlowerbd, CIupperbd, truelabel)
  vis_df <- as.data.frame(vis_df)
  names(vis_df) <- c('estimated', 'CIlowerbd', 'CIupperbd', 'truelabel')
  
   
  p1 <- ggplot(vis_df, aes(x = truelabel, y = estimated)) + 
		geom_point(size = 4) + geom_errorbar(aes(ymin = CIlowerbd, ymax = CIupperbd), size = 1) +
                xlab('Difference in Absolute Units') + ylab('Estimated Perceived Difference') +
		scale_color_discrete(name = 'Hypothesized\nFunction') #+
		#ggtitle('Discrete Estimation of Differences')

  p1

}
