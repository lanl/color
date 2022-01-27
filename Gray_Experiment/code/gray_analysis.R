## Investigation of Diminishing Returns in Gray Axis

rm(list = ls())
source('./code/functions/gray_functions.R')

wd <- './data/'

# Load in Data ------------------------------------------------------
d_form <- read.csv(paste(wd, 'gray_complete_data_release.csv', sep = ''))
d_form <- d_form %>%
    mutate(d1 = abs(Ls - Lt1),
           d2 = abs(Ls - Lt2))

# Significance Testing ----------------------------------------------
# Format the data
d_form_sig <- d_form %>% 
    group_by(Ls, Lt1, Lt2, d1, d2) %>%
    summarise(c_test2 = sum(R), tot = n()) %>%
    mutate(dbar = 0.5*(d1 + d2), dd = abs(d1 - d2))

d_form_sig$correct <- 0
d_form_sig <- d_form_sig %>%
    mutate(correct = replace(correct, d1 < d2, c_test2),
           correct = replace(correct, d1 > d2, tot - c_test2)) %>%
    group_by(dbar, dd) %>%
    summarise(acc = sum(correct) / sum(tot))

d_form_sig$dum_dd5 <- 0
d_form_sig[which(abs(d_form_sig$dd) == 5), 'dum_dd5'] <- 1

d_form_sig$dum_dd10 <- 0
d_form_sig[which(abs(d_form_sig$dd) == 10), 'dum_dd10'] <- 1

summary(d_form_sig$acc)

# Linear Regression
# Scaled
lm_out_sc <- lm(scale(acc) ~ dum_dd5 + dum_dd10 + scale(dbar), data = d_form_sig)
summary(lm_out_sc)

# Unscaled
lm_out <- lm(acc ~ dum_dd5 + dum_dd10 + dbar, data = d_form_sig)
summary(lm_out)

# Visualizing Accuracy
preg <- ggplot(d_form_sig, aes(x = dbar, y = acc, color = as.factor(abs(dd)))) + 
    geom_point() + geom_smooth(method = 'lm', se = TRUE) + theme_bw() +
    ylab('Accuracy') + xlab(expression('Average Difference,'~bar(d))) + 
    scale_color_discrete(name = expression(Delta*d))

preg

# MLE Analysis ------------------------------------------------------
# Set Parameters
my_sd <- 1
k <- 4

# Perform the Estimation
mle_out <- cv_boot_gray_gandf(d_form)
#save(mle_out, file = 'mle_out.Rdata')

# Visualize Accuracy
p_acc <- mle_out$master_accs %>% group_by(Ls, model) %>%
    summarise(GMacc = mean(Macc),
              SDacc = sd(Macc),
              SEacc = sd(Macc)/sqrt(n()),
              MEacc = 1.96*sd(Macc)/sqrt(n())) %>%
    ggplot(aes(x = Ls, y = GMacc, color = model, linetype = model)) + geom_point() + geom_line() +
        scale_color_manual(breaks = c('baseline', 'g only', 
                                      'g and f (spline)', 
                                      'g and f (polynomial)', 
                                      'g and f (log)', 
                                      'g and f (sine)'),
                          values = c('black', '#777777', '#1b9e77', '#d95f02', '#7570b3', '#e7298a'),
                          name = 'Model') +
        scale_linetype_manual(breaks = c('baseline', 'g only', 
                                      'g and f (spline)', 
                                      'g and f (polynomial)', 
                                      'g and f (log)', 
                                      'g and f (sine)'),
                              values = c('solid','dotted', 'dashed', 'longdash','dotdash','twodash'), name = 'Model')  +
        geom_errorbar(aes(ymin = GMacc - MEacc, ymax = GMacc + MEacc), width = 2) + theme_bw() +
        ylab('Average Accuracy for Predicting Responses') + xlab('L* Value of the Standard')

p_acc_sparse <- mle_out$master_accs %>% group_by(model) %>%
    summarise(GMacc = mean(Macc),
              SDacc = sd(Macc),
              SEacc = sd(Macc)/sqrt(n()),
              MEacc = 1.96*sd(Macc)/sqrt(n())) %>%
    ggplot(aes(x = factor(model, c('baseline', 'g only', 
                                      'g and f (spline)', 
                                      'g and f (polynomial)', 
                                      'g and f (log)', 
                                      'g and f (sine)')), y = GMacc)) + 
        geom_bar(stat = 'identity', fill = NA, color = 'black') +
        geom_errorbar(aes(ymin = GMacc - MEacc, ymax = GMacc + MEacc), width = 0.75) + theme_bw() +
        ylab('Average Accuracy for Predicting Responses') + xlab('Model') + 
        scale_y_continuous(limits = c(0.75, 0.925), oob = rescale_none)

# Visualizing g(x)
g_x_pars <- seq(0, 100, length.out = k+1)
g_x_plot <- seq(0, 100, length.out = 1000)

master_g_pars_only <- mle_out$master_g_pars %>% select(-cv_i, -bs_i)

master_g_line <- data.frame()

for(i in 1:nrow(master_g_pars_only)){
    y_tmp <- g_transform(g_x_pars, master_g_pars_only[i,], g_x_plot)
    tmp <- as.data.frame(cbind(g_x_plot, y_tmp))
    master_g_line <- rbind(master_g_line, tmp)
}

g_plot_df <- master_g_line %>% group_by(g_x_plot) %>%
    summarise(g_mean = mean(y_tmp),
              g_low = quantile(y_tmp, 0.025),
              g_upp = quantile(y_tmp, 0.975)) 

p_g <- g_plot_df %>% 
    ggplot(aes(x = g_x_plot, y = g_mean)) + 
        geom_ribbon(aes(ymin = g_low, ymax = g_upp), alpha = 0.2) + geom_line() +
        geom_segment(x = 0, y = 0, xend = 100, yend = max(g_plot_df$g_mean), linetype = 'dashed', size = .5) +
        ylab(expression('Perceived Strength,'~Psi)) + xlab('Absolute Strength, L*') + theme_bw()

# Using g(x) to transform data
agg_g_pars <- mle_out$master_g_pars %>% select(-cv_i, -bs_i) %>% colMeans()
# get the x pars for f spline
d_form_gtrans <- d_form %>% 
    mutate(PsiS = 100*g_transform(g_x_pars, agg_g_pars, Ls)/max(agg_g_pars),
           PsiT1 = 100*g_transform(g_x_pars, agg_g_pars, Lt1)/max(agg_g_pars),
           PsiT2 = 100*g_transform(g_x_pars, agg_g_pars, Lt2)/max(agg_g_pars),
           d1 = abs(PsiS - PsiT1),
           d2 = abs(PsiS - PsiT2)) 

mind <- min(c(d_form_gtrans$d1, d_form_gtrans$d2))
maxd <- max(c(d_form_gtrans$d1, d_form_gtrans$d2))

f_x_pars <- seq(mind, maxd, length.out = k+1)
f_x_plot <- seq(mind, maxd, length.out = 1000)

# f(x), spline
master_f_spl_pars_only <- mle_out$master_f_spl_pars %>% select(-cv_i, -bs_i)

master_f_spl_line <- data.frame()

for(i in 1:nrow(master_f_spl_pars_only)){
    y_tmp <- f_transform_spl(f_x_pars, master_f_spl_pars_only[i,], f_x_plot)
    tmp <- as.data.frame(cbind(f_x_plot, y_tmp))
    master_f_spl_line <- rbind(master_f_spl_line, tmp)
}

f_plot_spl_df <- master_f_spl_line %>% group_by(f_x_plot) %>%
    summarise(f_mean = mean(y_tmp),
              f_low = quantile(y_tmp, 0.025),
              f_upp = quantile(y_tmp, 0.975)) 

f_plot_spl_df$model <- 'spline'

# f(x), polynomial
fpar_low <- quantile(mle_out$master_f_ply_pars$fpar, 0.005)
fpar_upp <- quantile(mle_out$master_f_ply_pars$fpar, 0.995)

amp_low <- quantile(mle_out$master_f_ply_pars$amp, 0.005)
amp_upp <- quantile(mle_out$master_f_ply_pars$amp, 0.995)

master_f_ply_pars_only <- mle_out$master_f_ply_pars %>% select(-cv_i, -bs_i) %>%
    filter(fpar > fpar_low,
           fpar < fpar_upp,
           amp > amp_low,
           amp < amp_upp)

master_f_ply_line <- data.frame()

for(i in 1:nrow(master_f_ply_pars_only)){
    amp_i <- master_f_ply_pars_only[i, 'amp']
    fpar_i <- master_f_ply_pars_only[i, 'fpar']
    y_tmp <- f_transform_ply(amp_i, fpar_i, f_x_plot/maxd)
    
    tmp <- as.data.frame(cbind(f_x_plot, y_tmp))
    master_f_ply_line <- rbind(master_f_ply_line, tmp)
}

f_plot_ply_df <- master_f_ply_line %>% group_by(f_x_plot) %>%
    summarise(f_mean = mean(y_tmp),
              f_low = quantile(y_tmp, 0.025),
              f_upp = quantile(y_tmp, 0.975)) 

f_plot_ply_df$model <- 'polynomial'
f_plot_ply_df[1,1:4] <- 0

# f(x), log
fpar_low <- quantile(mle_out$master_f_log_pars$fpar, 0.005)
fpar_upp <- quantile(mle_out$master_f_log_pars$fpar, 0.995)

amp_low <- quantile(mle_out$master_f_log_pars$amp, 0.005)
amp_upp <- quantile(mle_out$master_f_ply_pars$amp, 0.995)

master_f_log_pars_only <- mle_out$master_f_log_pars %>% select(-cv_i, -bs_i) %>%
    filter(fpar > fpar_low,
           fpar < fpar_upp,
           amp > amp_low,
           amp < amp_upp)

master_f_log_line <- data.frame()

for(i in 1:nrow(master_f_log_pars_only)){
    amp_i <- master_f_log_pars_only[i, 'amp']
    fpar_i <- master_f_log_pars_only[i, 'fpar']
    y_tmp <- f_transform_log(amp_i, fpar_i, f_x_plot/maxd)
    
    tmp <- as.data.frame(cbind(f_x_plot, y_tmp))
    master_f_log_line <- rbind(master_f_log_line, tmp)
}

f_plot_log_df <- master_f_log_line %>% group_by(f_x_plot) %>%
    summarise(f_mean = mean(y_tmp),
              f_low = quantile(y_tmp, 0.025),
              f_upp = quantile(y_tmp, 0.975)) 

f_plot_log_df$model <- 'log'
f_plot_log_df[1,1:4] <- 0

# f(x), sine
fpar_low <- quantile(mle_out$master_f_sin_pars$fpar, 0.005)
fpar_upp <- quantile(mle_out$master_f_sin_pars$fpar, 0.995)

amp_low <- quantile(mle_out$master_f_sin_pars$amp, 0.005)
amp_upp <- quantile(mle_out$master_f_sin_pars$amp, 0.995)

master_f_sin_pars_only <- mle_out$master_f_sin_pars %>% select(-cv_i, -bs_i) %>%
    filter(fpar > fpar_low,
           fpar < fpar_upp,
           amp > amp_low,
           amp < amp_upp)

master_f_sin_line <- data.frame()

for(i in 1:nrow(master_f_sin_pars_only)){
    amp_i <- master_f_sin_pars_only[i, 'amp']
    fpar_i <- master_f_sin_pars_only[i, 'fpar']
    y_tmp <- f_transform_sin(amp_i, fpar_i, f_x_plot/maxd)
    
    tmp <- as.data.frame(cbind(f_x_plot, y_tmp))
    master_f_sin_line <- rbind(master_f_sin_line, tmp)
}

f_plot_sin_df <- master_f_sin_line %>% group_by(f_x_plot) %>%
    summarise(f_mean = mean(y_tmp),
              f_low = quantile(y_tmp, 0.025),
              f_upp = quantile(y_tmp, 0.975)) 

f_plot_sin_df$model <- 'sine'
f_plot_sin_df[1,1:4] <- 0

# Aggregating all f(x)
f_plot_df  <- f_plot_spl_df %>%
    bind_rows(f_plot_ply_df, f_plot_log_df, f_plot_sin_df)

# Visualizing f(x), spline
p_f_spline <- f_plot_df %>% filter(model == 'spline') %>%
ggplot(aes(x = f_x_plot, y = f_mean)) +
    geom_ribbon(aes(ymin = f_low, ymax = f_upp), fill = '#1b9e77', alpha = 0.1, linetype = 'blank') +
    geom_line(color = '#1b9e77') + theme_bw() + xlab(expression('Difference in Perceived Strengths,'~Delta*Psi)) +
    ylab('Perceived Difference') + geom_segment(x = 0, xend =  max(f_plot_spl_df$f_x_plot), 
                                                y = 0, yend =  max(f_plot_spl_df$f_mean), color = 'black',
                                               linetype = 'dashed', size = .5)

# Comparing other f(x) functions to spline
p_f_comp <- f_plot_df %>% filter(model != 'spline') %>%
ggplot(aes(x = f_x_plot, y = f_mean, color = model, fill = model)) +
    geom_ribbon(aes(ymin = f_low, ymax = f_upp), alpha = 0.1, linetype = 'blank') +
    geom_line() + theme_bw() + xlab(expression('Difference in Perceived Strengths,'~Delta*Psi)) +
    ylab('Perceived Difference') +
        scale_color_manual(breaks = c('polynomial', 
                                      'log', 
                                      'sine'),
                          values = c('#d95f02', '#7570b3', '#e7298a'),
                          name = 'Model')  + 
        scale_fill_manual(breaks = c('polynomial', 
                                      'log', 
                                      'sine'),
                          values = c('#d95f02', '#7570b3', '#e7298a'),
                          name = 'Model')  + 

geom_line(data = f_plot_spl_df, aes(x = f_x_plot, y = f_mean), color = 'black', inherit.aes = FALSE)
