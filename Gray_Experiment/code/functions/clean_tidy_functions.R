## Data Cleaning and Tidying Functions

library(dplyr)
library(reshape2)
library(stringr)

# First Reformatting ------------------------
first_form <- function(d_raw){
  d_raw <- d_raw[3:nrow(d_raw),]
  d_raw <- d_raw %>%
    filter(Finished == 'TRUE')
  d_raw
}


# Checking for Disclosed CVD ----------------
cvd_disclose <- function(d_raw){
  d_cvd_disc <- d_raw %>%
    filter(cvd == 'NO')
  d_cvd_disc
}


# Checking for CVD ------------------------
extract_correct <- function(string){
  if (grepl('null', string, fixed = T)) {
    0
  } else {
    c <- str_extract(string, pattern = '(\\d)+')
    as.numeric(c)
  }
}

cvd_check <- function(d_cvd_disc, cut = 2){
  d_cvd_test <- d_cvd_disc %>% #[,c(5:27,324)] 
    dplyr::select(c('ResponseId', starts_with('cvd'))) %>%
    dplyr::select(-cvd)
  
  d_cvd_test_melt <- melt(d_cvd_test, id.vars = 'ResponseId')
  d_cvd_test_melt <- d_cvd_test_melt %>%
    filter(value != '')

  d_cvd_test_melt$corr <- sapply(d_cvd_test_melt$variable,
                                 extract_correct)
  d_cvd_test_melt$point <- d_cvd_test_melt$value !=
    d_cvd_test_melt$corr

  d_cvd_test_res <- d_cvd_test_melt %>%
    group_by(ResponseId) %>%
    summarise(score = sum(as.numeric(point))) %>%
    filter(score < cut)

  d_cvd_clean <- d_cvd_test_res %>%
    left_join(d_cvd_disc, by = 'ResponseId')

  d_cvd_test_rejects <- d_cvd_test_melt %>%
    group_by(ResponseId) %>%
    summarise(score = sum(as.numeric(point))) %>%
    filter(score >= cut)

  my_list <- list('clean_data' = d_cvd_clean,
                  'cvd_reject' = d_cvd_test_rejects)

  my_list
}


# Formatting Data ----------------
formatting_data <- function(dict, d_cvd_clean, ntrial){
  names(dict)[4] <- 'value'
  dict_small <- dict[,c(6,7)]
  dict_small <- dict_small %>% filter(value != 'None')
  
  d_cvd_clean_form <- d_cvd_clean %>%
    dplyr::select(c(ResponseId, num_range("s", 0:ntrial),
                    num_range("c", 0:ntrial),
                    num_range("a", 0:ntrial),
                    num_range("b", 0:ntrial)))
  
  d_melt <- melt(d_cvd_clean_form,
                 id.vars = 'ResponseId')
  
  # D Choice
  d_choice <- d_melt %>%
    filter(grepl('c', variable))
  d_choice$trial <- str_split_fixed(d_choice$variable, 'c', 2)[,2]
  
  # D Trial Levels
  d_trial <- d_melt %>% 
    anti_join(d_choice, by = 'variable')
  
  d_trial <- d_trial %>%
    left_join(dict_small, by = 'value')
  
  d_trial$GrayPatch <- str_sub(d_trial$GrayPatch, 
                               5, str_length(d_trial$GrayPatch)-4)
  
  d_trial$trial <- str_split_fixed(d_trial$variable, '[sab]', 2)[,2]
  
  
  d_stand <- d_trial %>% filter(grepl('s', variable))
  d_stand <- d_stand %>% dplyr::select(ResponseId, trial, 'Ls' = GrayPatch)
  
  
  d_test1 <- d_trial %>% filter(grepl('a', variable))
  d_test1 <- d_test1 %>% dplyr::select(ResponseId, trial, 'Lt1' = GrayPatch)
  
  
  d_test2 <- d_trial %>% filter(grepl('b', variable))
  d_test2 <- d_test2 %>% dplyr::select(ResponseId, trial, 'Lt2' = GrayPatch)
  
  
  d_choice <- d_choice %>% dplyr::select(ResponseId, trial, 'R' = value)
  
  
  d_form <- d_stand %>%
    left_join(d_test1, by = c('ResponseId', 'trial')) %>%
    left_join(d_test2, by = c('ResponseId', 'trial')) %>%
    left_join(d_choice, by = c('ResponseId', 'trial')) #%>%
    #dplyr::select(ResponseId, Ls, Lt1, Lt2, R)
  
  d_form$Ls <- as.numeric(d_form$Ls)
  d_form$Lt1 <- as.numeric(d_form$Lt1)
  d_form$Lt2 <- as.numeric(d_form$Lt2)
  
  d_form$R <- as.numeric(as.factor(d_form$R))
  d_form$R <- d_form$R - 1 
  
  d_form
  
}
