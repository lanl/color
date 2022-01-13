## Hypothesis and Spline Analysis Functions
## Author: Emily Stark
## Date: 12-10-2020

# Full Notebook: hypothesis_spline_proposed_analysis.Rmd
# Qualtrics_full folder (Native on ES machine)

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(spacesXYZ)

## Homogeneity Check (plot only, options for bar or line) -------------------------------
homogen_check <- function(f1_data, type = 'bar'){ # Must have Ls, Lt1, Lt2, R
  
  if (type == 'line') {
    p1 <- f1_data %>%
      mutate(d1 = abs(Ls - Lt1),
             d2 = abs(Ls - Lt2),
             dd = d1 - d2,
             P = ifelse(dd < 0, R, -1*R+1)) %>%
      arrange(d1, d2) %>%
      mutate(triad = paste('d1_', d1, '_d2_', d2, sep = '')) %>%
      
      group_by(Ls, triad, d1 ,d2, dd) %>%
      summarise(choice2 = sum(P),
                total = n(),
                prop = sum(P)/n(),
                low = prop.test(sum(P),n())$conf.int[1],
                high = prop.test(sum(P),n())$conf.int[2]) %>%
      ungroup() %>%
      ggplot(aes(x = triad, y = prop, color = as.factor(Ls), 
                 group = as.factor(Ls))) +
        geom_point() + geom_line() +
        # geom_bar(stat = 'identity', position = 'dodge') +
        # geom_errorbar(aes(ymin = low, ymax = high), position = 'dodge') +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.75)) + 
        xlab('Triad') + ylab('Proportion Correct') + 
      scale_color_discrete(name = 'Center') +
        ggtitle('Check for Homogeneity') 
  }
  
  if (type == 'grid') {
    p1 <- f1_data %>%
      mutate(d1 = abs(Ls - Lt1),
             d2 = abs(Ls - Lt2),
             dd = d1 - d2,
             P = ifelse(dd < 0, R, -1*R+1)) %>%
      arrange(d1, d2) %>%
      mutate(triad = paste('d1_', d1, '_d2_', d2, sep = '')) %>%
      
      group_by(Ls, triad, d1 ,d2, dd) %>%
      summarise(choice2 = sum(P),
                total = n(),
                prop = sum(P)/n(),
                low = prop.test(sum(P),n())$conf.int[1],
                high = prop.test(sum(P),n())$conf.int[2]) %>%
      ungroup() %>%
      ggplot(aes(x = as.factor(Ls), y = prop, 
                 group = as.factor(Ls))) + facet_grid(d1 ~ dd, labeller = label_both) +
      
      geom_bar(stat = 'identity', position = 'dodge') +
      geom_errorbar(aes(ymin = low, ymax = high), position = 'dodge') +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.75)) + 
      xlab('Center') + ylab('Proportion Correct') + 
      scale_fill_discrete(name = 'Center') +
      ggtitle('Check for Homogeneity')
  }
  
  if (type == 'bar') {
    p1 <- f1_data %>%
      mutate(d1 = abs(Ls - Lt1),
             d2 = abs(Ls - Lt2),
             dd = d1 - d2,
             P = ifelse(dd < 0, R, -1*R+1)) %>%
      arrange(d1, d2) %>%
      mutate(triad = paste('d1_', d1, '_d2_', d2, sep = '')) %>%
      
     group_by(Ls, triad, d1 ,d2, dd) %>%
      summarise(choice2 = sum(P),
                total = n(),
                prop = sum(P)/n(),
                low = prop.test(sum(P),n())$conf.int[1],
                high = prop.test(sum(P),n())$conf.int[2]) %>%
      ungroup() %>%
      ggplot(aes(x = triad, y = prop, fill = as.factor(Ls), 
                 group = as.factor(Ls))) +
      # geom_point() + geom_line() +
      geom_bar(stat = 'identity', position = 'dodge') +
      geom_errorbar(aes(ymin = low, ymax = high), position = 'dodge') +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.75)) + 
      xlab('Triad') + ylab('Proportion Correct') + 
      scale_fill_discrete(name = 'Center') +
      ggtitle('Check for Homogeneity')
  }
  p1 
}



## Hypothesis Testing - Discovery -------------------------------------------------------
hyp_disc <- function(sim_data, homo = FALSE){
  if (homo == TRUE) {
    p1 <- sim_data %>%
      mutate(d1 = abs(Ls - Lt1),
             d2 = abs(Ls - Lt2),
             dbar = 0.5*(d1+d2),
             ddif = d1 - d2) %>%
      filter(dbar != 0) %>%
      filter(ddif != 0) %>%
      
      group_by(ddif, dbar) %>%
      summarise(totaln = n(),
                choos2 = sum(R)) %>%
      mutate(prop2 = choos2/totaln,
             prop1 = 1 - prop2,
             
             prop2 = replace(prop2, ddif < 0, 1 - prop2),
             prop1 = replace(prop1, ddif < 0, 1 - prop1),
             ddifa = abs(ddif)) %>%
      ungroup() %>%
      
      group_by(ddifa, dbar) %>%
      summarise(propc = mean(prop1)) %>%
      ungroup() %>%
      
      ggplot(aes(x = dbar, y = propc, group = ddifa, 
                 color = as.character(ddifa))) +
      geom_point(size = 3) + geom_smooth(method = 'lm', 
                                         formula = y~x, size = 1.5) +
      xlab('Average Difference in the Triad') +
      ylab('Proportion Correct Reponse') +
      scale_color_discrete(name = 'Delta D',
                           breaks = c('2.5', '5', '10'))
  }
  if (homo == FALSE) {
    p1 <- sim_data %>%
      mutate(d1 = abs(Ls - Lt1),
             d2 = abs(Ls - Lt2),
             dbar = 0.5*(d1+d2),
             ddif = d1 - d2) %>%
      filter(dbar != 0) %>%
      filter(ddif != 0) %>%
      
      group_by(ddif, dbar, Ls) %>%
      summarise(totaln = n(),
                choos2 = sum(R)) %>%
      mutate(prop2 = choos2/totaln,
             prop1 = 1 - prop2,
             
             prop2 = replace(prop2, ddif < 0, 1 - prop2),
             prop1 = replace(prop1, ddif < 0, 1 - prop1),
             ddifa = abs(ddif)) %>%
      ungroup %>%
      
      group_by(ddifa, dbar, Ls) %>%
      summarise(propc = mean(prop1)) %>%
      ungroup %>%
      
      ggplot(aes(x = dbar, y = propc, group = ddifa, 
                 color = as.character(ddifa))) +
      facet_wrap(vars(Ls), 
                 labeller = 
                   as_labeller(function(string, 
                                        suffix = "Center: ") paste0(
                                          suffix, string))) +
      geom_point(size = 3) + geom_smooth(method = 'lm', 
                                         formula = y~x, size = 1.5,
                                         se = FALSE) +
      xlab('Average Difference in the Triad') +
      ylab('Proportion Correct Reponse') +
      scale_color_discrete(name = 'Delta D',
                           breaks = c('2.5', '5', '10'))
  }
  p1
}


## Hypothesis Testing - Inference (Reg) -------------------------------------------------
hyp_reg <- function(sim_data, d.incl = F){
  if(d.incl){
  reg_df <- sim_data %>%
    mutate(dbar = 0.5*(d1+d2),
           ddif = d1 - d2) %>%
    filter(dbar != 0) %>%
    filter(ddif != 0) %>%
    
    group_by(ddif, dbar, Lt1) %>%
    summarise(totaln = n(),
              choos2 = sum(R)) %>%
    mutate(prop2 = choos2/totaln,
           prop1 = 1 - prop2,
           
           prop2 = replace(prop2, ddif < 0, 1 - prop2),
           prop1 = replace(prop1, ddif < 0, 1 - prop1),
           deltaD = abs(ddif)) %>%
    ungroup() %>%
    
    group_by(deltaD, dbar, Lt1) %>%
    summarise(propc = mean(prop1)) %>%
    ungroup() } else{
    reg_df <- sim_data %>%
    mutate(d1 = abs(Ls - Lt1),
           d2 = abs(Ls - Lt2),
           dbar = 0.5*(d1+d2),
           ddif = d1 - d2) %>%
    filter(dbar != 0) %>%
    filter(ddif != 0) %>%
    
    group_by(ddif, dbar, Lt1) %>%
    summarise(totaln = n(),
              choos2 = sum(R)) %>%
    mutate(prop2 = choos2/totaln,
           prop1 = 1 - prop2,
           
           prop2 = replace(prop2, ddif < 0, 1 - prop2),
           prop1 = replace(prop1, ddif < 0, 1 - prop1),
           deltaD = abs(ddif)) %>%
    ungroup() %>%
    
    group_by(deltaD, dbar, Lt1) %>%
    summarise(propc = mean(prop1)) %>%
    ungroup() 
    }
     
  reg_df$del5 <- ifelse(reg_df$deltaD == 5, 1, 0)
  reg_df$del10 <- ifelse(reg_df$deltaD == 10, 1, 0)
  
  my_out <- lm(propc ~ dbar + del5*dbar + del10*dbar +
                       del5 + del10, data = reg_df)

  #my_out <- lm(propc ~ dbar + deltaD + dbar*deltaD, data = reg_df)
  #my_out <- lm(propc ~ dbar, data = reg_df)
  summary(my_out)
  #table(reg_df$deltaD)
}

hyp_reg_nonhom <- function(sim_data){
  reg_df <- sim_data %>%
    mutate(d1 = abs(Ls - Lt1),
           d2 = abs(Ls - Lt2),
           dbar = 0.5*(d1+d2),
           ddif = d1 - d2) %>%
    filter(dbar != 0) %>%
    filter(ddif != 0) %>%
    
    group_by(ddif, dbar, Lt1, Ls) %>%
    summarise(totaln = n(),
              choos2 = sum(R)) %>%
    mutate(prop2 = choos2/totaln,
           prop1 = 1 - prop2,
           
           prop2 = replace(prop2, ddif < 0, 1 - prop2),
           prop1 = replace(prop1, ddif < 0, 1 - prop1),
           deltaD = abs(ddif)) %>%
    ungroup() %>%
    
    group_by(deltaD, dbar, Lt1, Ls) %>%
    summarise(propc = mean(prop1)) %>%
    ungroup() 
     
  reg_df$del5 <- ifelse(reg_df$deltaD == 5, 1, 0)
  reg_df$del10 <- ifelse(reg_df$deltaD == 10, 1, 0)
  
  my_out <- lm(propc ~ dbar + del5*dbar + del10*dbar +
                       del5 + del10 + Ls, data = reg_df)

  #my_out <- lm(propc ~ dbar + deltaD + dbar*deltaD, data = reg_df)
  #my_out <- lm(propc ~ dbar, data = reg_df)
  summary(my_out)
  #table(reg_df$deltaD)
}

hyp_reg_nonhom_s <- function(sim_data){
  reg_df <- sim_data %>%
    mutate(d1 = abs(Ls - Lt1),
           d2 = abs(Ls - Lt2),
           dbar = 0.5*(d1+d2),
           ddif = d1 - d2) %>%
    filter(dbar != 0) %>%
    filter(ddif != 0) %>%
    
    group_by(ddif, dbar, Lt1, Ls) %>%
    summarise(totaln = n(),
              choos2 = sum(R)) %>%
    mutate(prop2 = choos2/totaln,
           prop1 = 1 - prop2,
           
           prop2 = replace(prop2, ddif < 0, 1 - prop2),
           prop1 = replace(prop1, ddif < 0, 1 - prop1),
           deltaD = abs(ddif)) %>%
    ungroup() %>%
    
    group_by(deltaD, dbar, Lt1, Ls) %>%
    summarise(propc = mean(prop1)) %>%
    ungroup() 
     
  reg_df$del5 <- ifelse(reg_df$deltaD == 5, 1, 0)
  reg_df$del10 <- ifelse(reg_df$deltaD == 10, 1, 0)
  
  my_out <- lm(scale(propc) ~ scale(dbar) + del5*scale(dbar) + del10*scale(dbar) +
                       del5 + del10 + scale(Ls), data = reg_df)

  #my_out <- lm(propc ~ dbar + deltaD + dbar*deltaD, data = reg_df)
  #my_out <- lm(propc ~ dbar, data = reg_df)
  summary(my_out)
  #table(reg_df$deltaD)
}

## Hypothesis Testing - Inference (AOV) -------------------------------------------------
hyp_aov <- function(sim_data){
  aov_df <- sim_data %>%
    mutate(d1 = abs(Ls - Lt1),
           d2 = abs(Ls - Lt2),
           dbar = 0.5*(d1+d2),
           ddif = d1 - d2,
           center = as.factor(Ls)) %>%
    filter(dbar != 0) %>%
    filter(ddif != 0) %>%
    
    group_by(ddif, dbar, Lt1) %>%
    summarise(totaln = n(),
              choos2 = sum(R)) %>%
    mutate(prop2 = choos2/totaln,
           prop1 = 1 - prop2,
           
           prop2 = replace(prop2, ddif < 0, 1 - prop2),
           prop1 = replace(prop1, ddif < 0, 1 - prop1),
           deltaD = abs(ddif)) %>%
    ungroup %>%
    
    group_by(deltaD, dbar, Lt1) %>%
    summarise(propc = mean(prop1)) %>%
    ungroup %>%
    mutate(deltaD = as.factor(deltaD),
           dbar = as.factor(dbar),
           Lt1 = as.factor(Lt1))
    
  
  my_aov <- aov(propc ~ Lt1 + deltaD + dbar, 
                data = aov_df)
  
  summary(my_aov)
}

## Delta E Star 2000 Funtion ------------------------------------------------------------
DEstar2000 <- function(stand, test){
  cielab_stand <- cbind(stand, rep(0, length(stand)), rep(0, length(stand)))
  cielab_test <- cbind(test, rep(0, length(test)), rep(0, length(test)))
  DeltaE(cielab_stand, cielab_test, metric = '2000')
}

## Likelihood Function ------------------------------------------------------------------
continuous.lklh.freeMax <- function(data, par){
  my_par <- c(0, cumsum(abs(par)))
  f_pars <- my_par
  my_sd <- my_sd
    
  mind <- min(c(data$d1, data$d2))
  maxd <- max(c(data$d1, data$d2))
  
  x_pars <- seq(mind, maxd, length.out = length(f_pars))
  
  data %>%
    mutate(m1 = splinefun(x_pars, f_pars, 'monoH.FC')(d1), # transform d1
           m2 = splinefun(x_pars, f_pars, 'monoH.FC')(d2), # transform d2
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


## Multiple MLE Approximations ----------------------------------------------------------
multi_mle_better <- function(d_form_mle, k, my_sd){
    f_df <- data.frame()

    NLLs <- c()
    
    pb <- txtProgressBar(min = 0, max =100,
                         style = 3, width = 100, char= '=')
    for(i in 1:100){
        init_pars <- c(runif(k, 0, 1/(k)))
        test <- optim(par = init_pars,
                      fn = continuous.lklh.freeMax,
                      data = d_form_mle)

        tmp <- abs(test$par)
        tmp <- c(0, cumsum(tmp))

        f_df <- rbind(f_df, tmp)
        names(f_df) <- paste('V', seq(1,(k+1)), sep = '')

        NLLs <- c(NLLs, test$value)
        setTxtProgressBar(pb, i)
    }
    return(list('f_df' = f_df, 'NLLs' = NLLs))
}


