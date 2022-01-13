MLE_averaged <- function(MLE_output_df, L_vec){
  L_vec_min <- min(L_vec)
  L_vec_max <- max(L_vec)
  L_vec_ran <- L_vec_max - L_vec_min
  
  L_vec_cor <- -L_vec + L_vec_ran + 2*L_vec_min
  
  estimated_scaled <- sweep(MLE_output_df, MARGIN = 1, L_vec_cor, '*')
  estimated <- colSums(estimated_scaled) / sum(L_vec)

  return(estimated)
}


macadam_frac_exp <- function(par){
  if(par < mac_l) {par <- mac_l}
  if(par > mac_h) {par <- mac_h}
  m <- max(dis_df$dis_x)
  rmse_df <- dis_df %>% 
    mutate(dis_x = dis_x/m,
           dis_y_hat = dis_x^par,
           se = (dis_y_hat - dis_y)^2) %>%
    summarise(rmse = sqrt(mean(se))) 
  return(rmse_df[1,1])
}


helm_log <- function(par){
  m <- max(dis_df$dis_x)
  rmse_df <- dis_df %>% 
    mutate(dis_x = dis_x/m,
           dis_y_hat = log(par*dis_x + 1, 10)/log(par + 1, 10),
           dis_y_hat = replace(dis_y_hat, dis_y_hat == -Inf, 0),
           se = (dis_y_hat - dis_y)^2) %>%
    summarise(rmse = sqrt(mean(se))) 
  return(rmse_df[1,1])
}

izmailov_sin <- function(par){
  m <- max(dis_df$dis_x)
  rmse_df <- dis_df %>% 
    mutate(dis_x = dis_x/m,
           dis_y_hat = sin(dis_x*pi/(2*par))/sin(pi/(2*par)),

           se = (dis_y_hat - dis_y)^2) %>%
    summarise(rmse = sqrt(mean(se))) 
  return(rmse_df[1,1])
}

bujack_circle <- function(par){
  m <- max(dis_df$dis_x)
  rmse_df <- dis_df %>% 
    mutate(dis_x = dis_x/m,
           dis_y_hat = -par + 1 + sqrt(par^2 + 2*par*dis_x - dis_x^2 - 2*par + 1),
           se = (dis_y_hat - dis_y)^2) %>%
    summarise(rmse = sqrt(mean(se))) 
  return(rmse_df[1,1])
}
