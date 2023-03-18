## Analyzing In-Person Data
# Emily Teti
library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(reshape2)

set.seed(1234)

## Load in Data ------------------------------------------------------
setwd('C:/Users/353384/Projects/color')
data <- read.csv('inPerson_clean_data.csv')
names(data)[5] <- 'R'

data <- data %>% filter(Ls >= 50)

## Baseline, G, and F estimation functions ---------------------------
source('./functions/in_person_functions.R')


## Leave-One-Out CV to Establish Diminishing Returns -----------------
mle_out <- cv_loo_gray_gandf(data)

## Compare Accuracy --------------------------------------------------
p_acc <- mle_out$master_accs %>% group_by(model) %>%
  summarise(GMacc = mean(Macc),
            SDacc = sd(Macc),
            SEacc = sd(Macc)/sqrt(n()),
            MEacc = 1.96*sd(Macc)/sqrt(n()))
p_acc %>%
  ggplot(aes(x = factor(model, c('baseline', 'g only', 
                                 'g and f (spline)')), y = GMacc)) + 
  geom_bar(stat = 'identity', fill = NA, color = 'black') +
  geom_errorbar(aes(ymin = GMacc - MEacc, ymax = GMacc + MEacc), width = 0.75) + 
  theme_bw() +
  ylab('Average Accuracy for Predicting Responses') + xlab('Model') + 
  scale_y_continuous(limits = c(0.550, 0.675), oob = rescale_none)

lm_out <- lm (Macc~model, data = mle_out$master_accs)
aov_out <- aov(lm_out)
summary(aov_out)

TukeyHSD(aov_out)

## Visualize Perceptual Scale and Difference Scaling Function --------
g_x_pars <- seq(20, 100, length.out = k+1)
g_x_plot <- seq(20, 100, length.out = 1000)

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

g_plot_df %>% 
  ggplot(aes(x = g_x_plot, y = g_mean)) + 
  geom_ribbon(aes(ymin = g_low, ymax = g_upp), alpha = 0.2) + 
  geom_line() +
  geom_segment(x = 20, y = 0, xend = 100, yend = max(g_plot_df$g_mean), 
               linetype = 'dashed', size = .5) +
  ylab(expression('Perceived Strength,'~Psi)) + 
  xlab('Absolute Strength, L*') + theme_bw()

agg_g_pars <- mle_out$master_g_pars %>% select(-cv_i, -bs_i) %>% colMeans()

# get the x pars for f spline
d_form_gtrans <- data %>% 
  mutate(PsiS = 100*g_transform(g_x_pars, agg_g_pars, Ls)/max(agg_g_pars),
         PsiT1 = 100*g_transform(g_x_pars, agg_g_pars, Lt1)/max(agg_g_pars),
         PsiT2 = 100*g_transform(g_x_pars, agg_g_pars, Lt2)/max(agg_g_pars),
         d1 = abs(PsiS - PsiT1),
         d2 = abs(PsiS - PsiT2)) 

mind <- min(c(d_form_gtrans$d1, d_form_gtrans$d2))
maxd <- max(c(d_form_gtrans$d1, d_form_gtrans$d2))

f_x_pars <- seq(mind, maxd, length.out = k+1)
f_x_plot <- seq(mind, maxd, length.out = 1000)

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

f_plot_df  <- f_plot_spl_df

p_f_comp <- f_plot_df %>% 
  ggplot(aes(x = f_x_plot, y = f_mean)) +
  geom_ribbon(aes(ymin = f_low, ymax = f_upp), 
              alpha = 0.1, linetype = 'blank', fill = '#1b9e77') +
  geom_line(color = '#1b9e77') + theme_bw() + 
  xlab(expression('Difference in Perceived Strengths,'~Delta*Psi)) +
  ylab('Perceived Difference') #+
  # geom_line(data = f_plot_spl_df, 
  #           aes(x = f_x_plot, y = f_mean), color = 'black', inherit.aes = FALSE)

p_f_comp

## Comparing to the Crowdsourced Data ------------------------------------------
# Data and analysis come from github.com/lanl/color
data_cs <- read.csv('crowdsourced_clean_data.csv')
data_cs <- data_cs %>% filter(Ls >= 50)

## Comparing the Proportion Selecting the Lighter Test -------------------------
inPerson_summary <- data %>%
  mutate(d1 = abs(Ls - Lt1),
         d2 = abs(Ls - Lt2),
         dd = d1 - d2) %>%
  group_by(Ls, Lt1, Lt2, d1, d2, dd) %>%
  summarise(prop = mean(R))

crowdsourced_summary <- data_cs %>%
  mutate(d1 = abs(Ls - Lt1),
         d2 = abs(Ls - Lt2),
         dd = d1 - d2) %>%
  group_by(Ls, Lt1, Lt2, d1, d2, dd) %>%
  summarise(prop_cs = mean(R))

prop_comp <- inPerson_summary %>% left_join(crowdsourced_summary)
prop_comp %>% ggplot(aes(x = prop, y = prop_cs)) +
  geom_smooth(method = 'lm', color = 'black') + geom_point(aes(color = as.factor(dd))) + 
  scale_color_brewer(palette = 'PiYG', name = expression(paste(Delta, 'd'))) + 
  xlab('Proportion Selecting Lighter (in person)') +
  ylab('Proportion Selecting Lighter (crowdsourced)') +
  theme_bw()

cor.test(inPerson_summary$prop, crowdsourced_summary$prop_cs)


# Comparison of "Accuracy" -----------------------------------------------------
data_acc <- data
data_acc$correct <- 0
data_acc <- data_acc %>% mutate(Ls = as.numeric(Ls),
                        Lt1 = as.numeric(Lt1),
                        Lt2 = as.numeric(Lt2)) %>%
  mutate(d1 = abs(Ls - Lt1),
         d2 = abs(Ls - Lt2),
         dd = d1 - d2)


data_acc[which((data_acc$dd > 0) & (data_acc$R == 0)), 'correct'] <- 1
data_acc[which((data_acc$dd < 0) & (data_acc$R == 1)), 'correct'] <- 1

data_acc_cs <- data_cs
data_acc_cs$correct <- 0
data_acc_cs <- data_acc_cs %>% mutate(Ls = as.numeric(Ls),
                                Lt1 = as.numeric(Lt1),
                                Lt2 = as.numeric(Lt2)) %>%
  mutate(d1 = abs(Ls - Lt1),
         d2 = abs(Ls - Lt2),
         dd = d1 - d2)


data_acc_cs[which((data_acc_cs$dd > 0) & (data_acc_cs$R == 0)), 'correct'] <- 1
data_acc_cs[which((data_acc_cs$dd < 0) & (data_acc_cs$R == 1)), 'correct'] <- 1

mean(data_acc$correct)
mean(data_acc_cs$correct)

data_acc <- data_acc %>% group_by(Ls, Lt1, Lt2) %>%
  summarise(`In Person` = mean(correct))

data_acc_cs <- data_acc_cs %>% group_by(Ls, Lt1, Lt2) %>%
  summarise(Crowdsourced = mean(correct))

data_acc_comp <- data_acc %>% left_join(data_acc_cs)

data_acc_comp <- data_acc_comp %>% melt(id.vars = c('Ls', 'Lt1', 'Lt2'))

# p1 <- ggplot(data_acc_comp, aes(x = value, fill = variable)) + geom_histogram(alpha = 0.75) +
#   theme_bw() + 
#   scale_fill_brewer(palette = 'Dark2', name = 'Study') +
#   xlab('Accuracy') + ylab('Count') 
#   # scale_color_brewer(palette = 'Dark2')
# 
# ggsave('acc_comp_study.png')

## Analyze Crowdsourced Data ---------------------------------------------------
source('./functions/gray_functions.R')

my_sd <- 1
k <- 4

data_cs <- data_cs %>%
  mutate(d1 = abs(Ls - Lt1),
         d2 = abs(Ls - Lt2))

# Perform the Estimation
mle_out_cs <- cv_boot_gray_gandf(data_cs, bs_sampsize = 3200)

## Visualize Perceptual Scale and Difference Scaling Function (cs) --------
master_g_pars_only_cs <- mle_out_cs$master_g_pars %>% select(-cv_i, -bs_i)

master_g_line_cs <- data.frame()

for(i in 1:nrow(master_g_pars_only_cs)){
  y_tmp <- g_transform(g_x_pars, 
                       as.vector(master_g_pars_only_cs[i,]), 
                       g_x_plot)
  tmp <- as.data.frame(cbind(g_x_plot, y_tmp))
  master_g_line_cs <- rbind(master_g_line_cs, tmp)
  
}

g_plot_df_cs <- master_g_line_cs %>% group_by(g_x_plot) %>%
  summarise(g_mean = mean(y_tmp),
            g_low = quantile(y_tmp, 0.025),
            g_upp = quantile(y_tmp, 0.975)) 

g_plot_df$Method <- 'In Person'
g_plot_df$g_low <- g_plot_df$g_low / max(g_plot_df$g_mean)
g_plot_df$g_upp <- g_plot_df$g_upp / max(g_plot_df$g_mean)
g_plot_df$g_mean <- g_plot_df$g_mean / max(g_plot_df$g_mean)
  
  
g_plot_df_cs$Method <- 'Crowdsourced'
g_plot_df_cs$g_low <- g_plot_df_cs$g_low / max(g_plot_df_cs$g_mean)
g_plot_df_cs$g_upp <- g_plot_df_cs$g_upp / max(g_plot_df_cs$g_mean)
g_plot_df_cs$g_mean <- g_plot_df_cs$g_mean / max(g_plot_df_cs$g_mean)

g_plot_df_comp <- rbind(g_plot_df, g_plot_df_cs)


g_plot_df_comp %>% 
  ggplot(aes(x = g_x_plot, y = g_mean, color = Method, fill = Method)) + 
  geom_ribbon(aes(ymin = g_low, ymax = g_upp), alpha = 0.2) + 
  geom_line() +
  ylab(expression('Perceived Strength,'~Psi)) + 
  xlab('Absolute Strength, L*') + theme_bw() +
  scale_color_manual(breaks = c('Crowdsourced', 'In Person'),
                     values = c('#d95f02', '#7570b3')) +
  scale_fill_manual(breaks = c('Crowdsourced', 'In Person'),
                    values = c('#d95f02', '#7570b3'))

agg_g_pars <- mle_out$master_g_pars %>% select(-cv_i, -bs_i) %>% colMeans()

# get the x pars for f spline
d_form_gtrans_cs <- data_cs %>% 
  mutate(PsiS = 100*g_transform(g_x_pars, agg_g_pars, Ls)/max(agg_g_pars),
         PsiT1 = 100*g_transform(g_x_pars, agg_g_pars, Lt1)/max(agg_g_pars),
         PsiT2 = 100*g_transform(g_x_pars, agg_g_pars, Lt2)/max(agg_g_pars),
         d1 = abs(PsiS - PsiT1),
         d2 = abs(PsiS - PsiT2)) 

mind <- min(c(d_form_gtrans$d1, d_form_gtrans$d2))
maxd <- max(c(d_form_gtrans$d1, d_form_gtrans$d2))

f_x_pars <- seq(mind, maxd, length.out = k+1)
f_x_plot <- seq(mind, maxd, length.out = 1000)

master_f_spl_pars_only_cs <- mle_out_cs$master_f_spl_pars %>% select(-cv_i, -bs_i)

master_f_spl_line_cs <- data.frame()

for(i in 1:nrow(master_f_spl_pars_only_cs)){
  y_tmp <- f_transform_spl(f_x_pars, master_f_spl_pars_only_cs[i,], f_x_plot)
  tmp <- as.data.frame(cbind(f_x_plot, y_tmp))
  master_f_spl_line_cs <- rbind(master_f_spl_line_cs, tmp)
}

f_plot_spl_df_cs <- master_f_spl_line_cs %>% group_by(f_x_plot) %>%
  summarise(f_mean = mean(y_tmp),
            f_low = quantile(y_tmp, 0.025),
            f_upp = quantile(y_tmp, 0.975)) 

f_plot_spl_df_cs$model <- 'spline'

f_plot_df_cs  <- f_plot_spl_df_cs
f_plot_df_cs$Method <- 'Crowdsourced'
f_plot_df_cs$f_low <- f_plot_df_cs$f_low / max(f_plot_df_cs$f_mean)
f_plot_df_cs$f_upp <- f_plot_df_cs$f_upp / max(f_plot_df_cs$f_mean)
f_plot_df_cs$f_mean <- f_plot_df_cs$f_mean / max(f_plot_df_cs$f_mean)

f_plot_df$Method <- 'In Person'
f_plot_df$f_low <- f_plot_df$f_low / max(f_plot_df$f_mean)
f_plot_df$f_upp <- f_plot_df$f_upp / max(f_plot_df$f_mean)
f_plot_df$f_mean <- f_plot_df$f_mean / max(f_plot_df$f_mean)
f_plot_df$f_x_plot <- max(f_plot_df_cs$f_x_plot)*f_plot_df$f_x_plot / max(f_plot_df$f_x_plot)

f_plot_df_comp <- rbind(f_plot_df_cs, f_plot_df)


f_plot_df_comp %>% 
  ggplot(aes(x = f_x_plot, y = f_mean, color = Method, fill = Method)) + 
  geom_ribbon(aes(ymin = f_low, ymax = f_upp), alpha = 0.2, color = NA) + 
  geom_line() +
  xlab(expression('Difference in Perceived Strengths,'~Delta*Psi)) +
  ylab('Perceived Difference') + theme_bw() +
  scale_color_manual(breaks = c('Crowdsourced', 'In Person'),
                     values = c('#d95f02', '#7570b3')) +
  scale_fill_manual(breaks = c('Crowdsourced', 'In Person'),
                    values = c('#d95f02', '#7570b3'))
