# Load in Packages
library(dplyr, warn.conflicts = FALSE)
library(reshape2)
library(stringr)
library(ggplot2)

## g and f transformations ==============================================================
g_transform <- function(x_pars, g_pars, x){
    splinefun(x_pars, g_pars, 'monoH.FC')(x)
}

f_transform_spl <- function(x_pars, f_pars, x){
    splinefun(x_pars, f_pars, 'monoH.FC')(x)
}

f_transform_ply <- function(amp, fpar, x){
    amp*x^fpar
}

f_transform_log <- function(amp, fpar, x){
    amp*log(fpar*x + 1, 10)/log(fpar + 1, 10)
}

f_transform_sin <- function(amp, fpar, x){
    amp*sin(x*pi/(2*(fpar+0.5)))/sin(pi/(2*(fpar+0.5)))
}

## likelihood functions =================================================================
# scaling L axis (baseline)
# need Ls, Lt1, Lt2
# par is my_max
continuous.lklh.freeMax.L <- function(data, par){
  my_max <- abs(par)
  
  data %>%
    mutate(m1 = abs(Ls - Lt1)/my_max, # transform d1
           m2 = abs(Ls - Lt2)/my_max, # transform d2
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

# g transformation only
# need Ls, Lt1, Lt2
# par is 4 control points for g spline
continuous.lklh.freeMax.g <- function(data, par){
  my_par <- c(0, cumsum(abs(par)))
  g_pars <- my_par
  my_sd <- my_sd
    
  minL <- min(c(data$Ls, data$Lt1, data$Lt2))
  maxL <- max(c(data$Ls, data$Lt1, data$Lt2))
  
  x_pars <- seq(minL, maxL, length.out = length(g_pars))
  
  data %>%
    mutate(PsiS = g_transform(x_pars, g_pars, Ls), 
           PsiT1 = g_transform(x_pars, g_pars, Lt1),
           PsiT2 = g_transform(x_pars, g_pars, Lt2),
           m1 = abs(PsiS - PsiT1), # transform d1
           m2 = abs(PsiS - PsiT2), # transform d2
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

# f transformation, spline
# needs d1, d2
# par is 4 control points for f spline
continuous.lklh.freeMax.f.spl <- function(data, par){
    
    my_par <- c(0, cumsum(abs(par)))
    f_pars <- my_par
    
    mind <- min(c(data$d1, data$d2))
    maxd <- max(c(data$d1, data$d2))
  
    x_pars <- seq(mind, maxd, length.out = length(f_pars))
  
    data %>%
        mutate(m1 = f_transform_spl(x_pars, f_pars, d1),       # transform d1
               m2 = f_transform_spl(x_pars, f_pars, d2),       # transform d2
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

# f transformation, polynomial
# needs d1, d2
# par is amp and fpar for f polynomial
continuous.lklh.freeMax.f.ply <- function(data, par){
    
    amp <- abs(par[1])
    fpar <- abs(par[2])
    dmax <- max(c(data$d1, data$d2))
  
    data %>%
        mutate(m1 = f_transform_ply(amp, fpar, d1/dmax),       # transform d1
               m2 = f_transform_ply(amp, fpar, d2/dmax),       # transform d2
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

# f transformation, log
# needs d1, d2
# par is amp and fpar for f logarithmic
continuous.lklh.freeMax.f.log <- function(data, par){
    
    amp <- abs(par[1])
    fpar <- abs(par[2])
    dmax <- max(c(data$d1, data$d2))
  
    data %>%
        mutate(m1 = f_transform_log(amp, fpar, d1/dmax),       # transform d1
               m2 = f_transform_log(amp, fpar, d2/dmax),       # transform d2
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

# f transformation, sine
# needs d1, d2
# par is amp and fpar for f sine
continuous.lklh.freeMax.f.sin <- function(data, par){
    
    amp <- abs(par[1])
    fpar <- abs(abs(par[2]))
    dmax <- max(c(data$d1, data$d2))
  
    data %>%
        mutate(m1 = f_transform_sin(amp, fpar, d1/dmax),       # transform d1
               m2 = f_transform_sin(amp, fpar, d2/dmax),       # transform d2
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

## Accuracy Functions ===================================================================
# baseline
acc_baseline <- function(d_form_test, my_max, cv_i, bs_i, test_size_per_triad){
    tmp_acc <- d_form_test %>% group_by(Ls, Lt1, Lt2, d1, d2) %>%
    # Observed responses selecting test 2 in test set
    summarise(Rt2_obs = sum(R)) %>% 
    # Estimated responses selecting test 2 in test set
    mutate(m1 = abs(Ls - Lt1)/my_max, # transform d1
           m2 = abs(Ls - Lt2)/my_max, # transform d2
           dm = m1 - m2,
           Rt2_hat = pnorm(m2 - m1, mean = 0, 
                           sd = my_sd)*test_size_per_triad,
    # Calculate the accuracy, save with the bootstrap and cv indices
           acc = 1 - abs(Rt2_obs - Rt2_hat)/test_size_per_triad) %>% 
    ungroup() %>% group_by(Ls) %>%
    summarise(Macc = mean(acc)) %>%
    mutate(cv_ind = cv_i,
           bs_ind = bs_i,
           model = 'baseline')
    
    return(tmp_acc)
}

# g only
acc_g_only <- function(d_form_test, x_pars, g_pars, cv_i, bs_i, test_size_per_triad){
    tmp_acc <- d_form_test %>% group_by(Ls, Lt1, Lt2, d1, d2) %>%
    # Observed responses selecting test 2 in test set
    summarise(Rt2_obs = sum(R)) %>% 
    # Estimated responses selecting test 2 in test set
    mutate(PsiS = g_transform(x_pars, g_pars, Ls), 
           PsiT1 = g_transform(x_pars, g_pars, Lt1),
           PsiT2 = g_transform(x_pars, g_pars, Lt2),
           m1 = abs(PsiS - PsiT1), # transform d1
           m2 = abs(PsiS - PsiT2), # transform d2
           dm = m1 - m2,
           Rt2_hat = pnorm(m2 - m1, mean = 0, 
                           sd = my_sd)*test_size_per_triad,
    # Calculate the accuracy, save with the bootstrap and cv indices
           acc = 1 - abs(Rt2_obs - Rt2_hat)/test_size_per_triad) %>% 
    ungroup() %>% group_by(Ls) %>%
    summarise(Macc = mean(acc)) %>%
    mutate(cv_ind = cv_i,
           bs_ind = bs_i,
           model = 'g only')
    
    return(tmp_acc)
}

# g and f, spline
acc_f_spl <- function(d_form_test, x_pars, f_pars, cv_i, bs_i, test_size_per_triad){
    tmp_acc <- d_form_test %>% group_by(Ls, d1, d2) %>%
    # Observed responses selecting test 2 in test set
    summarise(Rt2_obs = sum(R)) %>% 
    # Estimated responses selecting test 2 in test set
    mutate(m1 = f_transform_spl(x_pars, f_pars, d1),       # transform d1
           m2 = f_transform_spl(x_pars, f_pars, d2),       # transform d2
           dm = m1 - m2,
           Rt2_hat = pnorm(m2 - m1, mean = 0, 
                           sd = my_sd)*test_size_per_triad,
    # Calculate the accuracy, save with the bootstrap and cv indices
           acc = 1 - abs(Rt2_obs - Rt2_hat)/test_size_per_triad) %>% 
    ungroup() %>% group_by(Ls) %>%
    summarise(Macc = mean(acc)) %>%
    mutate(cv_ind = cv_i,
           bs_ind = bs_i,
           model = 'g and f (spline)')
    
    return(tmp_acc)
}

# g and f, polynomial
acc_f_ply <- function(d_form_test, amp, fpar, cv_i, bs_i, test_size_per_triad){
    dmax <- max(c(d_form_test$d1, d_form_test$d2))

    tmp_acc <- d_form_test %>% group_by(Ls, d1, d2) %>%
    # Observed responses selecting test 2 in test set
    summarise(Rt2_obs = sum(R)) %>% 
    # Estimated responses selecting test 2 in test set
    mutate(m1 = f_transform_ply(amp, fpar, d1/dmax),       # transform d1
           m2 = f_transform_ply(amp, fpar, d2/dmax),       # transform d2
           dm = m1 - m2,
           Rt2_hat = pnorm(m2 - m1, mean = 0, 
                           sd = my_sd)*test_size_per_triad,
    # Calculate the accuracy, save with the bootstrap and cv indices
           acc = 1 - abs(Rt2_obs - Rt2_hat)/test_size_per_triad) %>% 
    ungroup() %>% group_by(Ls) %>%
    summarise(Macc = mean(acc)) %>%
    mutate(cv_ind = cv_i,
           bs_ind = bs_i,
           model = 'g and f (polynomial)')
    
    return(tmp_acc)
}

# g and f, log
acc_f_log <- function(d_form_test, amp, fpar, cv_i, bs_i, test_size_per_triad){
    dmax <- max(c(d_form_test$d1, d_form_test$d2))
    
    tmp_acc <- d_form_test %>% group_by(Ls, d1, d2) %>%
    # Observed responses selecting test 2 in test set
    summarise(Rt2_obs = sum(R)) %>% 
    # Estimated responses selecting test 2 in test set
    mutate(m1 = f_transform_log(amp, fpar, d1/dmax),       # transform d1
           m2 = f_transform_log(amp, fpar, d2/dmax),       # transform d2
           dm = m1 - m2,
           Rt2_hat = pnorm(m2 - m1, mean = 0, 
                           sd = my_sd)*test_size_per_triad,
    # Calculate the accuracy, save with the bootstrap and cv indices
           acc = 1 - abs(Rt2_obs - Rt2_hat)/test_size_per_triad) %>% 
    ungroup() %>% group_by(Ls) %>%
    summarise(Macc = mean(acc)) %>%
    mutate(cv_ind = cv_i,
           bs_ind = bs_i,
           model = 'g and f (log)')
    
    return(tmp_acc)
}

# g and f, sine
acc_f_sin <- function(d_form_test, amp, fpar, cv_i, bs_i, test_size_per_triad){
    dmax <- max(c(d_form_test$d1, d_form_test$d2))
    
    tmp_acc <- d_form_test %>% group_by(Ls, d1, d2) %>%
    # Observed responses selecting test 2 in test set
    summarise(Rt2_obs = sum(R)) %>% 
    # Estimated responses selecting test 2 in test set
    mutate(m1 = f_transform_sin(amp, fpar, d1/dmax),       # transform d1
           m2 = f_transform_sin(amp, fpar, d2/dmax),       # transform d2
           dm = m1 - m2,
           Rt2_hat = pnorm(m2 - m1, mean = 0, 
                           sd = my_sd)*test_size_per_triad,
    # Calculate the accuracy, save with the bootstrap and cv indices
           acc = 1 - abs(Rt2_obs - Rt2_hat)/test_size_per_triad) %>% 
    ungroup() %>% group_by(Ls) %>%
    summarise(Macc = mean(acc)) %>%
    mutate(cv_ind = cv_i,
           bs_ind = bs_i,
           model = 'g and f (sine)')
    
    return(tmp_acc)
}

## Cross Validation and Bootstrapping ===================================================
cv_boot_gray_gandf <- function(d_form, cv_N = 10, bs_N = 100,
                               test_size_per_triad = 50, bs_sampsize = 48000){
    options(dplyr.summarise.inform = FALSE)
    options(dplyr.join.inform = FALSE)

    
    # Create data frames to store results in
    master_L_par <- data.frame()
    master_g_pars <- data.frame()
    master_f_spl_pars <- data.frame()
    master_f_ply_pars <- data.frame()
    master_f_log_pars <- data.frame()
    master_f_sin_pars <- data.frame()
    
    master_accs <- data.frame()
    
    # Define the progress bar index
    pb_i <- 0

    pb <- txtProgressBar(min = 0, max = bs_n*cv_N,
                         style = 3, char= '=')
    
    # For the cross val index
    for(cv_i in 1:cv_N){
        # Test and train split
        d_form_test <- d_form %>% group_by(Ls, Lt1, Lt2) %>% 
            slice_sample(n = test_size_per_triad) 
        d_form_train <- anti_join(d_form, d_form_test, by = c("trial", "Ls", "Lt1", "Lt2", "R", "d1", "d2"))
        
        # Bootstrap a training set from the full training set
        for(bs_i in 1:bs_N){
            d_form_train_boot <- sample_n(d_form_train, bs_sampsize)
            
            # Complete the MLE for L* ===========================================================================================
            init_pars <- c(runif(1, 0, 100))
            test <- optim(par = init_pars,
                          method = 'Brent',
                          fn = continuous.lklh.freeMax.L,
                          data = d_form_train_boot,
                          lower = 0.00001, upper = 100)
            
            my_max <- abs(test$par)
            row <- c(my_max, cv_i, bs_i)
            master_L_par <- rbind(master_L_par, row)
            names(master_L_par) <- c('my_max', 'cv_i', 'bs_i')
            
            # Get Baseline Accuracy
            tmp_acc <- acc_baseline(d_form_test, my_max, cv_i, bs_i, test_size_per_triad)
            master_accs <- rbind(master_accs, tmp_acc) 
                        
            # Complete the MLE for g ============================================================================================
            init_pars <- c(runif(k, 0, 1/(k)))
            test <- optim(par = init_pars,
                  fn = continuous.lklh.freeMax.g,
                  data = d_form_train_boot)
            
            # Get the control points for g
            tmp <- abs(test$par)
            g_pars <- c(0, cumsum(tmp))
            
            row <- c(g_pars, cv_i, bs_i)
            master_g_pars <- rbind(master_g_pars, row)
            names(master_g_pars) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'cv_i', 'bs_i')
            
            minL <- min(c(d_form_train_boot$Ls, 
                          d_form_train_boot$Lt1, 
                          d_form_train_boot$Lt2))
            maxL <- max(c(d_form_train_boot$Ls, 
                          d_form_train_boot$Lt1, 
                          d_form_train_boot$Lt2))
            x_pars <- seq(minL, maxL, length.out = length(g_pars))
            
            
            # Compute Accuracy for g only broken up by center
            tmp_acc <- acc_g_only(d_form_test, x_pars, g_pars, cv_i, bs_i, test_size_per_triad)
            master_accs <- rbind(master_accs, tmp_acc) 
            
            # Transform L values in the data but maintain the scale
            d_form_train_boot_g <- d_form_train_boot %>%
                mutate(PsiS = g_transform(x_pars, g_pars, Ls), 
                       PsiT1 = g_transform(x_pars, g_pars, Lt1),
                       PsiT2 = g_transform(x_pars, g_pars, Lt2),
                
                       PsiS = maxL*PsiS/max(g_pars),
                       PsiT1 = maxL*PsiT1/max(g_pars),
                       PsiT2 = maxL*PsiT2/max(g_pars),
                       d1 = abs(PsiS - PsiT1), 
                       d2 = abs(PsiS - PsiT2))
            
            d_form_test_g <- d_form_test %>%
                mutate(PsiS = g_transform(x_pars, g_pars, Ls), 
                       PsiT1 = g_transform(x_pars, g_pars, Lt1),
                       PsiT2 = g_transform(x_pars, g_pars, Lt2),
                
                       PsiS = maxL*PsiS/max(g_pars),
                       PsiT1 = maxL*PsiT1/max(g_pars),
                       PsiT2 = maxL*PsiT2/max(g_pars),
                       d1 = abs(PsiS - PsiT1), 
                       d2 = abs(PsiS - PsiT2))
            
            # Complete the MLE for f, spline ====================================================================================
            init_pars <- c(runif(k, 0, 1/(k)))
            test <- optim(par = init_pars,
                          fn = continuous.lklh.freeMax.f.spl,
                          data = d_form_train_boot_g)

            # Get the control points for f
            tmp <- abs(test$par)
            f_pars <- c(0, cumsum(tmp))

            mind <- min(c(d_form_train_boot_g$d1, d_form_train_boot_g$d2))
            maxd <- max(c(d_form_train_boot_g$d1, d_form_train_boot_g$d2))
            x_pars <- seq(mind, maxd, length.out = length(f_pars))

            row <- c(f_pars, cv_i, bs_i)
            master_f_spl_pars <- rbind(master_f_spl_pars, row)
            names(master_f_spl_pars) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'cv_i', 'bs_i')
            
            # Calculate the accuracy, broken up by center, for f spline
            tmp_acc <- acc_f_spl(d_form_test_g, x_pars, f_pars, cv_i, bs_i, test_size_per_triad)
            master_accs <- rbind(master_accs, tmp_acc)
                        
            # Complete the MLE for f, polynomial ================================================================================
            init_pars <- c(runif(2, 0, 1))
            test <- optim(par = init_pars,
                          fn = continuous.lklh.freeMax.f.ply,
                          data = d_form_train_boot_g)

            # Get the parameters
            amp <- test$par[1]
            fpar <- test$par[2]

            row <- c(amp, fpar, cv_i, bs_i)
            master_f_ply_pars <- rbind(master_f_ply_pars, row)
            names(master_f_ply_pars) <- c('amp', 'fpar', 'cv_i', 'bs_i')

            # Calculate the accuracy, broken up by center, for f polynomial
            tmp_acc <- acc_f_ply(d_form_test_g, amp, fpar, cv_i, bs_i, test_size_per_triad)
            master_accs <- rbind(master_accs, tmp_acc)
            
            
            # Complete the MLE for f, log =======================================================================================
            init_pars <- c(runif(2, 0, 1))
            test <- optim(par = init_pars,
                          fn = continuous.lklh.freeMax.f.log,
                          data = d_form_train_boot_g)

            # Get the parameters
            amp <- test$par[1]
            fpar <- test$par[2]

            row <- c(amp, fpar, cv_i, bs_i)
            master_f_log_pars <- rbind(master_f_log_pars, row)
            names(master_f_log_pars) <- c('amp', 'fpar', 'cv_i', 'bs_i')

            # Calculate the accuracy, broken up by center, for f log
            tmp_acc <- acc_f_log(d_form_test_g, amp, fpar, cv_i, bs_i, test_size_per_triad)
            master_accs <- rbind(master_accs, tmp_acc)
            
            # Complete the MLE for f, sine ======================================================================================
            init_pars <- c(runif(2, 0, 1))
            test <- optim(par = init_pars,
                          fn = continuous.lklh.freeMax.f.sin,
                          data = d_form_train_boot_g)

            # Get the parameters
            amp <- test$par[1]
            fpar <- abs(test$par[2])

            row <- c(amp, fpar, cv_i, bs_i)
            master_f_sin_pars <- rbind(master_f_sin_pars, row)
            names(master_f_sin_pars) <- c('amp', 'fpar', 'cv_i', 'bs_i')

            # Calculate the accuracy, broken up by center, for f sine
            tmp_acc <- acc_f_sin(d_form_test_g, amp, fpar, cv_i, bs_i, test_size_per_triad)
            master_accs <- rbind(master_accs, tmp_acc)
            
            
            pb_i <- pb_i + 1
            setTxtProgressBar(pb, pb_i)

        }
    }
    
    my_list <- list('master_L_par' = master_L_par,
                    'master_g_pars' = master_g_pars,
                    'master_f_spl_pars' = master_f_spl_pars,
                    'master_f_ply_pars' = master_f_ply_pars,
                    'master_f_log_pars' = master_f_log_pars,
                    'master_f_sin_pars' = master_f_sin_pars,    
                    'master_accs' = master_accs) 
    return(my_list)
    
}

