
rm(list = ls())
if(Sys.info()[["user"]] == "s1437006"){
  setwd("/Users/s1437006/Documents/Edinburgh_PhD_documents/Projects/Trinidad_2019/Results/")
} else{
  setwd("/Users/cammybeyts/Documents/Edinburgh_PhD_documents/Projects/Trinidad_2019/Results/")
}

####load packages####
library(ggplot2)
library(dplyr) 
library(tidyverse)
library(brms)
library(rstanarm)
library(shinystan)
library(bayesplot)
library(kableExtra)
library(gridExtra)
library(emmeans)
library(tidybayes)

####ggplot themes #####

#ggplot theme for plotting figures
theme_catplots <- function() {
  theme_classic() + # use a white background
    theme(
      axis.title.x = element_text(face = "bold", size =20, colour = "black"),
      axis.title.y = element_text(face = "bold", size = 20, colour = "black"),
      axis.text.y = element_text(face = "bold", size = 14),
      axis.text.x = element_text(face = "bold", size = 14),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      strip.text.x = element_text(face = "bold", size=12),
      strip.background = element_rect(colour=FALSE, fill=FALSE),
      legend.position="none"
    )
}


#### functions for BRMS output #####

#the folllowing functions are used to obtain specific parameters from multivariate model

#create raw values for conspecific and heterospeific fixed effects
##create function to remove intercept, making treatment comparisons easier##

post_int_consp <- function(data_pd, params, type) {
  if (type == "M") {
    feNocomp <- data_pd[, paste("M", params[1], "int", sep = "_")]
    feConsp <- data_pd[, paste("M", params[1], "treatconsp", sep = "_")]
    feHetero <- data_pd[, paste("M", params[1], "treathetero", sep = "_")]
    feRaw <- data.frame(feNocomp + feConsp) 
  }
  else if (type == "S") {
    feNocomp <- data_pd[, paste("S", params[1], "int", sep = "_")]
    feConsp <- data_pd[, paste("S", params[1], "treatconsp", sep = "_")]
    feHetero <- data_pd[, paste("S", params[1], "treathetero", sep = "_")]
    feRaw <- data.frame(feNocomp + feConsp) 
  }
  else if (type == "H") {
    feNocomp <- data_pd[, paste("H", params[1], "int", sep = "_")]
    feConsp <- data_pd[, paste("H", params[1], "treatconsp", sep = "_")]
    feHetero <- data_pd[, paste("H", params[1], "treathetero", sep = "_")]
    feRaw <- data.frame(feNocomp + feConsp) 
  }
  return(feRaw)
}

post_int_hetero <- function(data_pd, params, type) {
  if (type == "M") {
    feNocomp <- data_pd[, paste("M", params[1], "int", sep = "_")]
    feConsp <- data_pd[, paste("M", params[1], "treatconsp", sep = "_")]
    feHetero <- data_pd[, paste("M", params[1], "treathetero", sep = "_")]
    feRaw <- data.frame(feNocomp + feHetero) 
  }
  else if (type == "S") {
    feNocomp <- data_pd[, paste("S", params[1], "int", sep = "_")]
    feConsp <- data_pd[, paste("S", params[1], "treatconsp", sep = "_")]
    feHetero <- data_pd[, paste("S", params[1], "treathetero", sep = "_")]
    feRaw <- data.frame(feNocomp + feHetero) 
  }
  else if (type == "H") {
    feNocomp <- data_pd[, paste("H", params[1], "int", sep = "_")]
    feConsp <- data_pd[, paste("H", params[1], "treatconsp", sep = "_")]
    feHetero <- data_pd[, paste("H", params[1], "treathetero", sep = "_")]
    feRaw <- data.frame(feNocomp + feHetero) 
  }
  return(feRaw)
}


##create function to perform difference calculations

#no competition - conspecific treatment
post_diff_nocomp_consp <- function(data_pd, params, type) {
  if (type == "V") {
    vcNocomp <- data_pd[, paste("V_id", params[1], "Nocomp", sep = "_")]
    vcConsp <- data_pd[, paste("V_id", params[1], "Consp", sep = "_")]
    vcHetero <- data_pd[, paste("V_id", params[1], "Hetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcConsp))
  }
  else if (type == "M") {
    vcNocomp <- data_pd[, paste("M", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("M", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("M", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcConsp))
  }
  else if (type == "S") {
    vcNocomp <- data_pd[, paste("S", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("S", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("S", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcConsp))
  }
  else if (type == "H") {
    vcNocomp <- data_pd[, paste("H", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("H", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("H", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcConsp))
  }
  else if (type == "corr") {
    vcNocomp <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Nocomp",
      sep = "_"
    )]
    vcConsp <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Consp",
      sep = "_"
    )]
    vcHetero <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Hetero",
      sep = "_"
    )]
    vcDiff <- data.frame(vcNocomp - vcConsp)
  }
  return(vcDiff)
}

post_diff_nocomp_consp_R <- function(data_pd, params, type) {
  if (type == "R") {
    vcNocomp <- data_pd[, paste("R", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("R", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("R", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcConsp))
  }}

#no competition - heterospecific treatment
post_diff_nocomp_hetero <- function(data_pd, params, type) {
  if (type == "V") {
    vcNocomp <- data_pd[, paste("V_id", params[1], "Nocomp", sep = "_")]
    vcConsp <- data_pd[, paste("V_id", params[1], "Consp", sep = "_")]
    vcHetero <- data_pd[, paste("V_id", params[1], "Hetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcHetero))
  }
  else if (type == "M") {
    vcNocomp <- data_pd[, paste("M", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("M", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("M", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcHetero))
  }
  else if (type == "S") {
    vcNocomp <- data_pd[, paste("S", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("S", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("S", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcHetero))
  }
  else if (type == "H") {
    vcNocomp <- data_pd[, paste("H", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("H", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("H", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcHetero))
  }
  else if (type == "corr") {
    vcNocomp <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Nocomp",
      sep = "_"
    )]
    vcConsp <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Consp",
      sep = "_"
    )]
    vcHetero <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Hetero",
      sep = "_"
    )]
    vcDiff <- data.frame(vcNocomp - vcHetero)
  }
  return(vcDiff)
}

post_diff_nocomp_hetero_R <- function(data_pd, params, type) {
  if (type == "R") {
    vcNocomp <- data_pd[, paste("R", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("R", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("R", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcNocomp - vcHetero))
  }}

#conspecific - heterospecific treatment
post_diff_consp_hetero <- function(data_pd, params, type) {
  if (type == "V") {
    vcNocomp <- data_pd[, paste("V_id", params[1], "Nocomp", sep = "_")]
    vcConsp <- data_pd[, paste("V_id", params[1], "Consp", sep = "_")]
    vcHetero <- data_pd[, paste("V_id", params[1], "Hetero", sep = "_")]
    vcDiff <- abs(data.frame(vcConsp - vcHetero))
  }
  if (type == "M") {
    vcNocomp <- data_pd[, paste("M", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("M", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("M", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcConsp - vcHetero))
  }
  if (type == "S") {
    vcNocomp <- data_pd[, paste("S", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("S", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("S", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcConsp - vcHetero))
  }
  if (type == "H") {
    vcNocomp <- data_pd[, paste("H", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("H", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("H", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcConsp - vcHetero))
  }
  else if (type == "corr") {
    vcNocomp <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Nocomp",
      sep = "_"
    )]
    vcConsp <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Consp",
      sep = "_"
    )]
    vcHetero <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Hetero",
      sep = "_"
    )]
    vcDiff <- data.frame(vcConsp - vcHetero)
  }
  return(vcDiff)
}

post_diff_consp_hetero_R <- function(data_pd, params, type) {
  if (type == "R") {
    vcNocomp <- data_pd[, paste("R", params[1], "treatnocomp", sep = "_")]
    vcConsp <- data_pd[, paste("R", params[1], "treatconsp", sep = "_")]
    vcHetero <- data_pd[, paste("R", params[1], "treathetero", sep = "_")]
    vcDiff <- abs(data.frame(vcConsp - vcHetero))
  }}
####Get data into R####
read.data <- read.delim("Trinidad_2019_ACT_EXP_PRED.txt", header = TRUE)

#copy treatment column and rename
read.data <- read.data %>%
  mutate(NameTreatment=Treatment)

#convert order values 
read.data <- read.data %>%
  mutate(Treatment=recode(Treatment, NoCompetition = 0, Conspecific = 1, Heterospecific= 2)) %>% #treatments converted into 1, 2 and 3
  mutate(Rep=recode(Rep, 0, 1, 2, 3, 4, 5)) %>%  #reps 123456 recoded to 012345 %>%
  mutate(Scale_SVL = scale(SVL))

#set up data types
full.data1 <- read.data %>%
  mutate(Treatment=as.factor(Treatment)) %>%
  mutate(Set=as.factor(Set)) %>% 
  mutate(RepCon=as.numeric(Rep))

#for activity data: change values of 0 into the lowest value recorded for lognormal distribution model 
full.data1$ACT.Dist.Pixels[full.data1$ACT.Dist.Pixels %in% 0] <- 1.41

####brms chains####
n_thin <- 10
burnin <- 1000
pd_size <- 3000
n_chains <- 4
n_iter <- round(pd_size / n_chains) * n_thin + burnin


####brms - body size univariate ####

SVL_brms_2022a <- brm(
    SVL ~ Treatment + RepCon + (1|FocalNest) + (1 | TadpoleID),
  data = full.data1,
  family = gaussian(),
  chains = n_chains, iter = n_iter, warmup = burnin, thin = n_thin
)
#save(SVL_brms_2022a, file = "SVL_brms_2022a.rda")
#load(file = "SVL_brms_2022a.rda")
summary(SVL_brms_2022a)

#plot model outputs
SVL_brms_2022a_areas <- mcmc_areas(SVL_brms_2022a, 
                                regex_pars = c("^b_"), prob = 0.95, # 80% intervals
                                prob_outer = 0.95, # 95%
                                point_est = "mean")
SVL_brms_2022a_areas

#### get svl univariate results #####

#fixed effects
fixef_mod_svl <- fixef(SVL_brms_2022a, summary = FALSE)
#variance for tadpole ID
v_id_svl <- (VarCorr(SVL_brms_2022a, summary = FALSE)$TadpoleID$sd)^2
#variances for eggmass ID
v_eggmass_svl <- (VarCorr(SVL_brms_2022a, summary = FALSE)$FocalNest$sd)^2

column_names_svl <- c("Mean", 
                      "2.5%", 
                      "97.5%"
)

row_names_svl <- c("Intercept",
                   "Conspecific",
                   "Heterospecific",
                   "Trial",
                   "Tadpole ID",
                   "Egg Mass ID"
)

#get summary of fixed effects
fixef_mod_SVL <- fixef_mod_svl[,c(1:4)]
fixef_mod_SVL <- round(data.frame(colMeans(fixef_mod_SVL), HPDinterval(as.mcmc(fixef_mod_SVL))), digits = 3) 
colnames(fixef_mod_SVL) <- column_names_svl

#get summary of eggmass
v_eggmass_SVL <- v_eggmass_svl[,1]
v_eggmass_SVL <- round(data.frame(mean(v_eggmass_SVL), HPDinterval(as.mcmc(v_eggmass_SVL))), digits = 3)
colnames(v_eggmass_SVL) <- column_names_svl

#get summary of variances 
v_id_SVL <- v_id_svl
v_id_SVL <- round(data.frame(colMeans(v_id_SVL), HPDinterval(as.mcmc(v_id_SVL))), digits = 3) 
colnames(v_id_SVL) <- column_names_svl

#combine fixed effects, EggMassID and ID variance into one data frame
SVL_df <- rbind(fixef_mod_SVL, v_id_SVL, v_eggmass_SVL )
rownames(SVL_df) <- row_names_svl

#kable table of body size results
SVL_kable <- SVL_df  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "Body Size" = 3))
SVL_kable


#### SVL plot preparation ####

#create parameter names of treatments
params_names_svl <- c(
  "M_m_int",
  "M_m_treatconsp", 
  "M_m_treathetero"
)

##subset fixed effects from my model
#all parameters
m_svl_all <- as.data.frame(SVL_brms_2022a)
#just select the fixed effects from the list of parameters
m_svl <- m_svl_all[c(1:3)] 

colnames(m_svl) <- params_names_svl

#currently mean body size estimates of tadpoles in the heterospecific and conspecific treatments are compared to the intercept
#apply post_int_consp and post_int_hetero functions to obtain posterior estimates for the mean body size of tadpoles in the heterospecific and conspecific treatments.  The intercept will be the posterior estimate of the no competition treatment

m_svl <- m_svl %>%
  mutate(
    M_m_treatconsp = unlist(post_int_consp(., params = "m", type = "M")), #posterior estimate of overall scaled mean body size, conspecific treatment
    M_m_treathetero = unlist(post_int_hetero(., params = "m", type = "M")), #posterior estimate of overall scaled mean body size, heterospefific treatment
    M_m_treatnocomp = M_m_int #posterior estimate of overall scaled mean body size, no competition treatment
  )

#use the post_diff_nocomp_consp, post_diff_nocomp_hetero and post_diff_consp_hetero functions to obtain a mean estimate and 95% CI of the difference body size between treatments. 

##average body size difference  between treatments
#no competition - conspecific treatment
m_svl_m_nocomp_consp <- m_svl %>%
  mutate(
    d_M_fe_am_ncc = unlist(post_diff_nocomp_consp(., params = "m", type = "M"))
  )

#no competition - heterospecific treatment
m_svl_m_nocomp_hetero <- m_svl %>%
  mutate(
    d_M_fe_am_nch = unlist(post_diff_nocomp_hetero(., params = "m", type = "M"))
  )

#conspecific - heterospecific treatment
m_svl_m_consp_hetero <- m_svl %>%
  mutate(
    d_M_fe_am_ch = unlist(post_diff_consp_hetero(., params = "m", type = "M"))
  )

##combine data frames together to make one data frame for:
#difference in mean behaviour between treatments
m_svl_m_nocomp_consp_hetero  <- cbind(m_svl,
                                      dplyr::select(m_svl_m_nocomp_consp, d_M_fe_am_ncc),
                                      dplyr::select(m_svl_m_nocomp_hetero, d_M_fe_am_nch),
                                      dplyr::select(m_svl_m_consp_hetero, d_M_fe_am_ch)
                                      )

#### fixed effect plots - svl #####

#create a plot of the average scaled body size of tadpoles in each treatment 
#and the difference in average scaled body size of tadpoles between treatments. 

plot_SVL_feMeanDiff_nocomp_consp_hetero <- 
  ggplot(
    (mcmc_intervals_data(m_svl_m_nocomp_consp_hetero,
                         regex_pars = 
                           "^M_m_int|^M_m_treatconsp|^M_m_treathetero|^d_M_fe_am_ncc|^d_M_fe_am_nch|d_M_fe_am_ch",
                         prob = 0.8, # 80% intervals
                         prob_outer = 0.95, # 95%
                         point_est = "mean"
    ) %>% mutate(
      trait = rep(c(" "),6),
      treat = c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff")
    )
    ),
    aes(y = treat, color = treat)
  ) +
  geom_point(aes(x = m), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = ll, xmax = hh)) +
  geom_linerange(aes(xmin = l, xmax = h), lwd = 2) +
  labs(
    title = bquote("Mean body size and \n difference in body size between treatments"),
    x = "\n Posterior Estimate of body size (mm) \n with 80% (thin) 95% (thick) Credible Intervals",
    y = "\n Treatment"
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0.001, 0), linetype = "dotted") +
  theme_catplots() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
plot_SVL_feMeanDiff_nocomp_consp_hetero
#ggsave(plot_SVL_feMeanDiff_nocomp_consp_hetero, filename = "SVL_pot.png", height = 8, width = 12)




##### mulitvariate BRMS model ####

#activity assay model
bf_ACT <- bf(
  ACT.Dist.Pixels ~ Treatment + Scale_SVL + RepCon + (1|FocalNest) + (1 |a| gr(TadpoleID, by = Treatment))) +
  lf(sigma ~ Treatment) + lognormal()

#exploration assay model 
bf_EXP <- bf(
  EXP.Dist.Pixels ~ Treatment + Scale_SVL + RepCon + (1|FocalNest) + (1 |a| gr(TadpoleID, by = Treatment))) +
  lf(hu ~ Treatment, sigma ~ Treatment) + hurdle_lognormal()

#predation-risk assay model
bf_PRED <- bf(
  PRED.Dist.Pixels ~ Treatment + Scale_SVL + RepCon + (1|FocalNest) + (1 |a| gr(TadpoleID, by = Treatment))) +
  lf(hu ~ Treatment, sigma ~ Treatment) + hurdle_lognormal()

#combine activity, exploration and predation-risk assay models into one multivariate model
ACT_EXP_PRED_brms_2022a <- brm(
  bf_ACT + bf_EXP + bf_PRED + set_rescor(FALSE), 
  chains = n_chains, iter = n_iter, warmup = burnin, thin = n_thin,
  data = full.data1
)
summary(ACT_EXP_PRED_brms_2022a)
#save(ACT_EXP_PRED_brms_2022a, file = "ACT_EXP_PRED_brms_2022a.rda")
#load(file = "ACT_EXP_PRED_brms_2022a.rda")

#plot model outputs
ACT_EXP_PRED_brms_2022a_areas <- mcmc_areas(ACT_EXP_PRED_brms_2022a, 
                                         regex_pars = c("^sd_"), prob = 0.80, #use "^sd_" to look at mcmc_areas of variances, "^b_" for fixed effects "^cor_" for correlations
                                         prob_outer = 0.95, # 95%
                                         point_est = "mean")
ACT_EXP_PRED_brms_2022a_areas


#####get multivariate model results - var and fixed effects####

#fixed effects
fixef_mod <- fixef(ACT_EXP_PRED_brms_2022a, summary = FALSE)
#variances for tadpole ID
v_id <- (VarCorr(ACT_EXP_PRED_brms_2022a, summary = FALSE)$TadpoleID$sd)^2
#correalations
cor_id_list <- VarCorr(ACT_EXP_PRED_brms_2022a, summary = TRUE)$TadpoleID$cor
#variances for eggmass ID
v_eggmass <- (VarCorr(ACT_EXP_PRED_brms_2022a, summary = FALSE)$FocalNest$sd)^2

column_names <- c("Mean", 
                  "2.5%", 
                  "97.5%"
)


row_names_ACTdummy <- c("Hurdle_Intercept",
                        "Hurdle_Conspecific",
                        "Hurdle_Heterospecific"
)


row_names_rep <- c("No Competition", "Conspecific", "Heterospecific")

#get summary of fixed effects
fixef_mod_ACT <- fixef_mod[,c(1:2, 9:14)]
fixef_mod_ACT <- round(data.frame(colMeans(fixef_mod_ACT), HPDinterval(as.mcmc(fixef_mod_ACT))), digits = 3) 
fixef_mod_EXP <- fixef_mod[,c(3:5, 15:22)]
fixef_mod_EXP <- round(data.frame(colMeans(fixef_mod_EXP), HPDinterval(as.mcmc(fixef_mod_EXP))), digits = 3) 
fixef_mod_PRED <- fixef_mod[,c(6:8, 23:30)]
fixef_mod_PRED <- round(data.frame(colMeans(fixef_mod_PRED), HPDinterval(as.mcmc(fixef_mod_PRED))), digits = 3) 
colnames(fixef_mod_ACT) <- column_names
colnames(fixef_mod_EXP) <- column_names
colnames(fixef_mod_PRED) <- column_names

#get summary of eggmass
v_eggmass_ACT <- v_eggmass[,1]
v_eggmass_ACT <- round(data.frame(mean(v_eggmass_ACT), HPDinterval(as.mcmc(v_eggmass_ACT))), digits = 3)
v_eggmass_EXP <- v_eggmass[,2]
v_eggmass_EXP <- round(data.frame(mean(v_eggmass_EXP), HPDinterval(as.mcmc(v_eggmass_EXP))), digits = 3)
v_eggmass_PRED <- v_eggmass[,3]
v_eggmass_PRED <- round(data.frame(mean(v_eggmass_PRED), HPDinterval(as.mcmc(v_eggmass_PRED))), digits = 3)
colnames(v_eggmass_ACT) <- column_names
colnames(v_eggmass_EXP) <- column_names
colnames(v_eggmass_PRED) <- column_names

#get summary of variances 
v_id_ACT <- v_id[,c(1, 4, 7)]
v_id_ACT <- round(data.frame(colMeans(v_id_ACT), HPDinterval(as.mcmc(v_id_ACT, prob=0.95))), digits = 3) 
v_id_EXP <- v_id[,c(2, 5, 8)]
v_id_EXP <- round(data.frame(colMeans(v_id_EXP), HPDinterval(as.mcmc(v_id_EXP))), digits = 3) 
v_id_PRED <- v_id[,c(3, 6, 9)]
v_id_PRED <- round(data.frame(colMeans(v_id_PRED), HPDinterval(as.mcmc(v_id_PRED))), digits = 3) 
colnames(v_id_ACT) <- column_names
colnames(v_id_EXP) <- column_names
colnames(v_id_PRED) <- column_names

#combine results of fixed effects and variances into one dataframe for each assay
EXP_df <- rbind(fixef_mod_EXP, v_eggmass_EXP, v_id_EXP)
PRED_df <- rbind(fixef_mod_PRED, v_eggmass_PRED, v_id_PRED)
ACT_df <- rbind(fixef_mod_ACT, v_eggmass_ACT, v_id_ACT)

#Ensure paramaters of each assay are displayed in the same order
EXP_df <- EXP_df %>%
  slice(c(1, 4, 5, 6, 7, 12, 13, 14, 15, 2, 8, 9, 3, 10, 11)) 
PRED_df <- PRED_df %>%
  slice(c(1, 4, 5, 6, 7, 12, 13, 14, 15, 2, 8, 9, 3, 10, 11)) 
ACT_df <- ACT_df %>%
  slice(c(1, 3, 4, 5, 6, 9, 10, 11, 12, 2, 7, 8)) 

#activity assay does not have hurdel part of model, NA need to be added into the dataframe for these
ACT_dum <- data.frame(matrix(NA, nrow = 3, ncol = 3,
                             dimnames = list(NULL, paste0("ColumnName_", 1:3))) )
(ACT_dum) <- column_names
rownames(ACT_dum) <- row_names_ACTdummy
ACT_df <- rbind(ACT_df,ACT_dum)

#combine activity, exploration and predation dataframes together and label the row names. 
ACT_EXP_PRED_df <- cbind(ACT_df, EXP_df, PRED_df)
rownames(ACT_EXP_PRED_df) <- c("Intercept",
                               "Conspecific",
                               "Heterospecific",
                               "SVL",
                               "Trial",
                               "Egg Mass ID",
                               "VI No Competition",
                               "VI Conspecific",
                               "VI Heterospecific",
                               "Intercept ", #sigma intercept
                               "Conspecific ",#sigma consp
                               "Heterospecific ", #sigma hetero
                               "Intercept  ", #hurdle intercept
                               "Conspecific  ", #hurdle consp
                               "Heterospecific  " #hurdle hetero
)

#create kable table of the final dataframe
ACT_EXP_PRED_kable <- ACT_EXP_PRED_df  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "Activity" = 3, "Exploration" = 3, "Predation" = 3)) %>%
  pack_rows("Mean Model", 1, 9) %>%
  pack_rows("Residuals", 10, 12) %>%
  pack_rows("Hurdle Model", 13, 15)
ACT_EXP_PRED_kable

##### get multivariate model results - correaltions #####

#pairwise correalations between assays
cor_id_list <- VarCorr(ACT_EXP_PRED_brms_2022a, summary = FALSE)$TadpoleID$cor
cor_id_list_nocomp <- cor_id_list[, 1:3, 1:3] %>%
  as.data.frame.array()
cor_id_list_consp <- cor_id_list[, 4:6, 4:6] %>%
  as.data.frame.array()
cor_id_list_hetero <- cor_id_list[, 7:9, 7:9] %>%
  as.data.frame.array()

#pairewise correlations between assays, treatment differences
cor_id_list_nocomp_consp <- cor_id_list_nocomp - cor_id_list_consp
cor_id_list_nocomp_hetero <- cor_id_list_nocomp - cor_id_list_hetero
cor_id_list_consp_hetero <- cor_id_list_consp - cor_id_list_hetero

row_names_corr <-  c(
  paste(c("ActVI", "ExpVI", "PredVI"), rep("ActVI", 3), sep = "_"),
  paste(c("ActVI", "ExpVI", "PredVI"), rep("ExpVI", 3), sep = "_"),
  paste(c("ActVI", "ExpVI", "PredVI"), rep("PredVI", 3), sep = "_")
)

#obtain mean estimate and 95% CI of pairwise correalations between assays for each treatment
#no competitions treatment
cor_id_list_nocomp <-  round(data.frame(colMeans(cor_id_list_nocomp), HPDinterval(as.mcmc(cor_id_list_nocomp))), digits = 3) %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))
#conspeific treatment
cor_id_list_consp  <-  round(data.frame(colMeans(cor_id_list_consp), HPDinterval(as.mcmc(cor_id_list_consp))), digits = 3)  %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))
#heterospeific treatment
cor_id_list_hetero  <-  round(data.frame(colMeans(cor_id_list_hetero), HPDinterval(as.mcmc(cor_id_list_hetero))), digits = 3)  %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))

#obtain mean estimate and 95% CI of pairwise correalations between assays - treatment comparisons
#no competition - conspecific treatment
cor_id_list_nocomp_consp <-  round(data.frame(colMeans(cor_id_list_nocomp_consp), HPDinterval(as.mcmc(cor_id_list_nocomp_consp))), digits = 3) %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))
#no competition - heterospeific treatment
cor_id_list_nocomp_hetero  <-  round(data.frame(colMeans(cor_id_list_nocomp_hetero), HPDinterval(as.mcmc(cor_id_list_nocomp_hetero))), digits = 3)  %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))
#conspecific - heterospeific treatment
cor_id_list_consp_hetero  <-  round(data.frame(colMeans(cor_id_list_consp_hetero), HPDinterval(as.mcmc(cor_id_list_consp_hetero))), digits = 3)  %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))

#create single dataframe of pairwise correalations for all treatments
cor_nocomp_consp_hetero <- cbind(cor_id_list_nocomp, cor_id_list_consp, cor_id_list_hetero) 
#create single dataframe of pairwise correalations for all treatment comparisons
cor_nocomp_consp_nocomp_hetero_consp_hetero <- cbind(cor_id_list_nocomp_consp, cor_id_list_nocomp_hetero, cor_id_list_consp_hetero)

#create table of pairwise correalations for all treatments
cor_nocomp_consp_hetero_kable <- cor_nocomp_consp_hetero  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "No Competition" = 3, "Conspecific" = 3, "Heterospecific" = 3))
cor_nocomp_consp_hetero_kable

#create table of pairwise correalations for all treatment comparisons
cor_nocomp_consp_nocomp_hetero_consp_hetero_kable <- cor_nocomp_consp_nocomp_hetero_consp_hetero  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = T) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "No Competition-Conspecific" = 3, "No Competition-Heterospecific" = 3, "Conspecific-Heterospecific" = 3))
cor_nocomp_consp_nocomp_hetero_consp_hetero_kable


#### plot preparation: data subset, functions and dataframes for plots ####

#create parameter names
params_names <- c(
  "M_act_m_int",
  "S_act_v_int", 
  "M_exp_m_int",
  "S_exp_v_int", 
  "H_exp_h_int", 
  "M_pred_m_int", 
  "S_pred_v_int", 
  "H_pred_h_int", 
  "M_act_m_treatconsp", "M_act_m_treathetero", "act_m_SVL", "act_m_rep",
  "S_act_v_treatconsp", "S_act_v_treathetero",
  "M_exp_m_treatconsp", "M_exp_m_treathetero", "exp_m_SVL", "exp_m_rep",
  "S_exp_v_treatconsp", "S_exp_v_treathetero",
  "H_exp_h_treatconsp", "H_exp_h_treathetero",
  "M_pred_m_treatconsp", "M_pred_m_treathetero", "pred_m_SVL", "pred_m_rep",
  "S_pred_v_treatconsp", "S_pred_v_treathetero",
  "H_pred_h_treatconsp", "H_pred_h_treathetero",
  "V_n_am",
  "V_id_am_Nocomp",
  "V_id_em_Nocomp",
  "V_id_pm_Nocomp",
  "V_id_am_Consp",
  "V_id_em_Consp",
  "V_id_pm_Consp",
  "V_id_am_Hetero",
  "V_id_em_Hetero",
  "V_id_pm_Hetero",
  "V_n_em",
  "V_n_pm",
  
  paste("cor_id", rep("am", 2), c("em", "pm"), "Nocomp", sep = "_"),
  paste("cor_id", rep("em"), c("pm"), "Nocomp", sep = "_"),
  
  paste("cor_id", rep("am", 2), c("em", "pm"), "Consp", sep = "_"),
  paste("cor_id", rep("em"), c("pm"), "Consp", sep = "_"),
  
  
  paste("cor_id", rep("am", 2), c("em", "pm"), "Hetero", sep = "_"),
  paste("cor_id", rep("em"), c("pm"), "Hetero", sep = "_")
)

##subset data from model
mv_act_exp_pred_all <- as.data.frame(ACT_EXP_PRED_brms_2022a)
colnames(mv_act_exp_pred_all)
names(mv_act_exp_pred_all[c(1:51)])

#fixed effects
mv_act_exp_pred_fx <- mv_act_exp_pred_all[c(1:30)] 
#variance of tadpole ID and Nest ID
mv_act_exp_pred_var <- mv_act_exp_pred_all[,31:42]^2
#combine all parameters into one dataframe
mv_act_exp_pred <- cbind(mv_act_exp_pred_fx, mv_act_exp_pred_var, mv_act_exp_pred_corr)
colnames(mv_act_exp_pred) <- params_names

##use post_int_consp and post_int_hetero functions to obtain mean estimates of activity, exploration and predation behaviour - current model displays these estimates as comparisons to the interept##

mv_act_exp_pred <- mv_act_exp_pred %>%
  mutate(
    M_act_m_treatconsp = unlist(post_int_consp(., params = "act_m", type = "M")), #posterior estimate of overall mean activity, conspecific treatment
    M_exp_m_treatconsp = unlist(post_int_consp(., params = "exp_m", type = "M")),  #posterior estimate of overall mean exploration, conspecific treatment
    M_pred_m_treatconsp = unlist(post_int_consp(., params = "pred_m", type = "M")), #posterior estimate of overall mean predation, conspecific treatment
    S_act_v_treatconsp = unlist(post_int_consp(., params = "act_v", type = "S")), #posterior estimate of variance within individuals in activity, conspecific treatment
    S_exp_v_treatconsp = unlist(post_int_consp(., params = "exp_v", type = "S")), #posterior estimate of variance within individuals in exploration, conspecific treatment
    S_pred_v_treatconsp = unlist(post_int_consp(., params = "pred_v", type = "S")), #posterior estimate of variance within individuals in predation, conspecific treatment
    H_exp_h_treatconsp = unlist(post_int_consp(., params = "exp_h", type = "H")), #posterior estimate of hurdle model in exploration, conspecific treatment
    H_pred_h_treatconsp = unlist(post_int_consp(., params = "pred_h", type = "H")), #posterior estimate of hurdle model in predation, conspecific treatment
    M_act_m_treathetero = unlist(post_int_hetero(., params = "act_m", type = "M")), #posterior estimate of overall mean activity, heterospecific treatment
    M_exp_m_treathetero = unlist(post_int_hetero(., params = "exp_m", type = "M")), #posterior estimate of overall mean exploration, heterospecific treatment
    M_pred_m_treathetero = unlist(post_int_hetero(., params = "pred_m", type = "M")),#posterior estimate of overall mean predation, heterospecific treatment
    S_act_v_treathetero = unlist(post_int_hetero(., params = "act_v", type = "S")), #posterior estimate of variance within individuals in activity, heterospecific treatment
    S_exp_v_treathetero = unlist(post_int_hetero(., params = "exp_v", type = "S")), #posterior estimate of variance within individuals in exploration, heterospecific treatment
    S_pred_v_treathetero = unlist(post_int_hetero(., params = "pred_v", type = "S")), #posterior estimate of variance within individuals in predation, heterospecific treatment
    H_exp_h_treathetero = unlist(post_int_hetero(., params = "exp_h", type = "H")), #posterior estimate of hurdle model in exploration, heterospecific treatment
    H_pred_h_treathetero = unlist(post_int_hetero(., params = "pred_h", type = "H")), #posterior estimate of hurdle model in predation, heterospecific treatment
    M_act_m_treatnocomp = M_act_m_int, #posterior estimate of overall mean activity, no competition treatment
    M_exp_m_treatnocomp = M_exp_m_int, #posterior estimate of overall mean exploration, no competition treatment
    M_pred_m_treatnocomp = M_pred_m_int, #posterior estimate of overall mean predation, no competition treatment
    S_act_v_treatnocomp = S_act_v_int, #posterior estimate of variance within individuals in activity, no competition treatment
    S_exp_v_treatnocomp = S_exp_v_int, #posterior estimate of variance within individuals in exploration, no competition treatment
    S_pred_v_treatnocomp = S_pred_v_int,  #posterior estimate of variance within individuals in predation, no competition treatment
    H_exp_h_treatnocomp = H_exp_h_int, #posterior estimate of hurdle model in exploration, no competition treatment
    H_pred_h_treatnocomp = H_pred_h_int, #posterior estimate of hurdle model in predation, no competition treatment
  )

#calculate repeatabilties
rep_act_exp_pred <- mv_act_exp_pred %>%
  mutate(
    R_act_m_treatnocomp = mv_act_exp_pred$V_id_am_Nocomp/(mv_act_exp_pred$V_id_am_Nocomp + mv_act_exp_pred$S_act_v_treatnocomp), #repeatability of activity behaviour, no competition treatment 
    R_act_m_treatconsp = mv_act_exp_pred$V_id_am_Consp/(mv_act_exp_pred$V_id_am_Consp + mv_act_exp_pred$S_act_v_treatconsp), #repeatability of activity behaviour, conspecific treatment 
    R_act_m_treathetero = mv_act_exp_pred$V_id_am_Hetero/(mv_act_exp_pred$V_id_am_Hetero + mv_act_exp_pred$S_act_v_treathetero), #repeatability of activity behaviour, heterospecific treatment 
    R_exp_m_treatnocomp = mv_act_exp_pred$V_id_em_Nocomp/(mv_act_exp_pred$V_id_em_Nocomp + mv_act_exp_pred$S_exp_v_treatnocomp), #repeatability of exploration behaviour, no competition treatment 
    R_exp_m_treatconsp = mv_act_exp_pred$V_id_em_Consp/(mv_act_exp_pred$V_id_em_Consp + mv_act_exp_pred$S_exp_v_treatconsp), #repeatability of exploration behaviour, conspeific treatment 
    R_exp_m_treathetero = mv_act_exp_pred$V_id_em_Hetero/(mv_act_exp_pred$V_id_em_Hetero + mv_act_exp_pred$S_exp_v_treathetero), #repeatability of exploration behaviour, heterospecific treatment 
    R_pred_m_treatnocomp = mv_act_exp_pred$V_id_pm_Nocomp/(mv_act_exp_pred$V_id_pm_Nocomp + mv_act_exp_pred$S_pred_v_treatnocomp), #repeatability of predation risk-taking behaviour, no competition treatment 
    R_pred_m_treatconsp = mv_act_exp_pred$V_id_pm_Consp/(mv_act_exp_pred$V_id_pm_Consp + mv_act_exp_pred$S_pred_v_treatconsp),  #repeatability of predation risk-taking behaviour, conspecific treatment 
    R_pred_m_treathetero = mv_act_exp_pred$V_id_pm_Hetero/(mv_act_exp_pred$V_id_pm_Hetero + mv_act_exp_pred$S_pred_v_treathetero) #repeatability of predation risk-taking behaviour, heterospecific treatment 
  )

###Use post_diff functions to calculate difference in parameter estimates between treatments

#no competition - conspecific treatment
mv_act_exp_pred_m_nocomp_consp <- mv_act_exp_pred %>%
  mutate(
    d_V_id_am_ncc = unlist(post_diff_nocomp_consp(., params = "am", type = "V")), #difference in posterior estimate of overall activity behaviour 
    d_V_id_em_ncc = unlist(post_diff_nocomp_consp(., params = "em", type = "V")), #difference in posterior estimate of overall exploration behaviour 
    d_V_id_pm_ncc = unlist(post_diff_nocomp_consp(., params = "pm", type = "V")), #difference in posterior estimate of overall predation risk taking behaviour 
    d_M_fe_am_ncc = unlist(post_diff_nocomp_consp(., params = "act_m", type = "M")), #difference in posterior estimate of among individual variance in activity behaviour 
    d_M_fe_em_ncc = unlist(post_diff_nocomp_consp(., params = "exp_m", type = "M")), #difference in posterior estimate of among individual variance in exploration behaviour
    d_M_fe_pm_ncc = unlist(post_diff_nocomp_consp(., params = "pred_m", type = "M")), #difference in posterior estimate of among individual variance in predation risk taking behaviour
    d_S_fe_av_ncc = unlist(post_diff_nocomp_consp(., params = "act_v", type = "S")), #difference in posterior estimate of within individual variance in activity behaviour
    d_S_fe_ev_ncc = unlist(post_diff_nocomp_consp(., params = "exp_v", type = "S")), #difference in posterior estimate of within individual variance in exploration behaviour
    d_S_fe_pv_ncc = unlist(post_diff_nocomp_consp(., params = "pred_v", type = "S")), #difference in posterior estimate of within individual variance in predation risk taking behaviour
    d_H_fe_eh_ncc = unlist(post_diff_nocomp_consp(., params = "exp_h", type = "H")), #difference in posterior estimate of hurdle model exploration behaviour
    d_H_fe_ph_ncc = unlist(post_diff_nocomp_consp(., params = "pred_h", type = "H")) #difference in posterior estimate of hurdle model predation risk taking behaviour
  )


mv_act_exp_pred_m_nocomp_consp_R <- rep_act_exp_pred %>%
  mutate(
    d_R_id_am_ncc = unlist(post_diff_nocomp_consp_R(., params = "act_m", type = "R")), #difference in posterior estimate of the repeatability of activity behaviour 
    d_R_id_em_ncc = unlist(post_diff_nocomp_consp_R(., params = "exp_m", type = "R")), #difference in posterior estimate of the repeatability of exploration behaviour 
    d_R_id_pm_ncc = unlist(post_diff_nocomp_consp_R(., params = "pred_m", type = "R")) #difference in posterior estimate of the repeatability of predatory risk taking behaviour 
  )

#no competition - heterospecific treatment
mv_act_exp_pred_m_nocomp_hetero <- mv_act_exp_pred %>%
  mutate(
    d_V_id_am_nch = unlist(post_diff_nocomp_hetero(., params = "am", type = "V")),
    d_V_id_em_nch = unlist(post_diff_nocomp_hetero(., params = "em", type = "V")),
    d_V_id_pm_nch = unlist(post_diff_nocomp_hetero(., params = "pm", type = "V")),
    d_M_fe_am_nch = unlist(post_diff_nocomp_hetero(., params = "act_m", type = "M")),
    d_M_fe_em_nch = unlist(post_diff_nocomp_hetero(., params = "exp_m", type = "M")),
    d_M_fe_pm_nch = unlist(post_diff_nocomp_hetero(., params = "pred_m", type = "M")),
    d_S_fe_av_nch = unlist(post_diff_nocomp_hetero(., params = "act_v", type = "S")),
    d_S_fe_ev_nch = unlist(post_diff_nocomp_hetero(., params = "exp_v", type = "S")),
    d_S_fe_pv_nch = unlist(post_diff_nocomp_hetero(., params = "pred_v", type = "S")),
    d_H_fe_eh_nch = unlist(post_diff_nocomp_hetero(., params = "exp_h", type = "H")),
    d_H_fe_ph_nch = unlist(post_diff_nocomp_hetero(., params = "pred_h", type = "H")),
    d_cor_id_am_em_nch = unlist(post_diff_nocomp_hetero(., params = c("am", "em"), type = "corr")),
    d_cor_id_am_pm_nch = unlist(post_diff_nocomp_hetero(., params = c("am", "pm"), type = "corr")),
    d_cor_id_em_pm_nch = unlist(post_diff_nocomp_hetero(., params = c("em", "pm"), type = "corr"))
  )
#repeatability
mv_act_exp_pred_m_nocomp_hetero_R <- rep_act_exp_pred %>%
  mutate(
    d_R_id_am_nch = unlist(post_diff_nocomp_hetero_R(., params = "act_m", type = "R")),
    d_R_id_em_nch = unlist(post_diff_nocomp_hetero_R(., params = "exp_m", type = "R")),
    d_R_id_pm_nch = unlist(post_diff_nocomp_hetero_R(., params = "pred_m", type = "R"))
  )

#conspecific - heterospecific treatment
mv_act_exp_pred_m_consp_hetero <- mv_act_exp_pred %>%
  mutate(
    d_V_id_am_ch = unlist(post_diff_consp_hetero(., params = "am", type = "V")),
    d_V_id_em_ch = unlist(post_diff_consp_hetero(., params = "em", type = "V")),
    d_V_id_pm_ch = unlist(post_diff_consp_hetero(., params = "pm", type = "V")),
    d_M_fe_am_ch = unlist(post_diff_consp_hetero(., params = "act_m", type = "M")),
    d_M_fe_em_ch = unlist(post_diff_consp_hetero(., params = "exp_m", type = "M")),
    d_M_fe_pm_ch = unlist(post_diff_consp_hetero(., params = "pred_m", type = "M")),
    d_S_fe_av_ch = unlist(post_diff_consp_hetero(., params = "act_v", type = "S")),
    d_S_fe_ev_ch = unlist(post_diff_consp_hetero(., params = "exp_v", type = "S")),
    d_S_fe_pv_ch = unlist(post_diff_consp_hetero(., params = "pred_v", type = "S")),
    d_H_fe_eh_ch = unlist(post_diff_consp_hetero(., params = "exp_h", type = "H")),
    d_H_fe_ph_ch = unlist(post_diff_consp_hetero(., params = "pred_h", type = "H")),
    d_cor_id_am_em_ch = unlist(post_diff_consp_hetero(., params = c("am", "em"), type = "corr")),
    d_cor_id_am_pm_ch = unlist(post_diff_consp_hetero(., params = c("am", "pm"), type = "corr")),
    d_cor_id_em_pm_ch = unlist(post_diff_consp_hetero(., params = c("em", "pm"), type = "corr"))
  )
#repeatability
mv_act_exp_pred_m_consp_hetero_R <- rep_act_exp_pred %>%
  mutate(
    d_R_id_am_ch = unlist(post_diff_consp_hetero_R(., params = "act_m", type = "R")),
    d_R_id_em_ch = unlist(post_diff_consp_hetero_R(., params = "exp_m", type = "R")),
    d_R_id_pm_ch = unlist(post_diff_consp_hetero_R(., params = "pred_m", type = "R"))
  )

##combine data frames together to make one data frame for:
  #difference in overall mean, variance among individuals, variance within individuals, hurdle model between treatments
mv_act_exp_pred_m_nocomp_consp_hetero  <- cbind(mv_act_exp_pred, 
                                                dplyr::select(mv_act_exp_pred_m_nocomp_consp, d_V_id_am_ncc:d_cor_id_em_pm_ncc), 
                                                dplyr::select(mv_act_exp_pred_m_nocomp_hetero, d_V_id_am_nch:d_cor_id_em_pm_nch),
                                                dplyr::select(mv_act_exp_pred_m_consp_hetero, d_V_id_am_ch:d_cor_id_em_pm_ch)
                                                )
  #difference in repeatability between treatments
mv_act_exp_pred_s_nocomp_consp_hetero_R <- cbind(rep_act_exp_pred,
                                                 dplyr::select(mv_act_exp_pred_m_nocomp_consp_R, d_R_id_am_ncc:d_R_id_pm_ncc),
                                                 dplyr::select(mv_act_exp_pred_m_nocomp_hetero_R, d_R_id_am_nch:d_R_id_pm_nch),
                                                 dplyr::select(mv_act_exp_pred_m_consp_hetero_R, d_R_id_am_ch:d_R_id_pm_ch)
                                                 )  

#### variance among and within and repeatability plots ######

#create plot of variance among individuals in activity, exloratin and predation assays 
#plot include estimate of the among individual level variance in activity, exloratin and predation risk taking between treatments

plot_ACTEXPPRED_VarMeanDiff_nocomp_consp_hetero <- 
  ggplot(
    (mcmc_intervals_data(mv_act_exp_pred_m_nocomp_consp_hetero,
                         regex_pars = 
                           "^V_id_am_Nocomp|^V_id_am_Consp|V_id_am_Hetero|d_V_id_am_ncc|d_V_id_am_nch|d_V_id_am_ch|V_id_em_Nocomp|V_id_em_Consp|V_id_em_Hetero|d_V_id_em_ncc|d_V_id_em_nch|d_V_id_em_ch|V_id_pm_Nocomp|V_id_pm_Consp|V_id_pm_Hetero|d_V_id_pm_ncc|d_V_id_pm_nch|d_V_id_pm_ch",
                         prob = 0.8, # 80% intervals
                         prob_outer = 0.95, # 95%
                         point_est = "median"
    ) %>% mutate(
      trait = rep(c("Activity", "Exporation", "Predation"), 6),
      treat = rep(c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"), each = 3)
    )
    ),
    aes(y = treat, color = treat)
  ) +
  geom_point(aes(x = m), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = ll, xmax = hh)) +
  geom_linerange(aes(xmin = l, xmax = h), lwd = 2) +
  labs(
    title = bquote("Variance Among Individuals"),
    x = " ", # " \n Posterior Estimate with 80% (thin) 95% (thick) Credible Intervals",
    y = " " #"Treatment"
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=10))
plot_ACTEXPPRED_VarMeanDiff_nocomp_consp_hetero

#create plot of variance within individuals in activity, exloratin and predation assays 
#plot include estimate of the within individual level variance in activity, exloratin and predation risk taking between treatments

plot_ACTEXPPRED_feSigmaDiff_nocomp_consp_hetero <- 
  ggplot(
    (mcmc_intervals_data(mv_act_exp_pred_m_nocomp_consp_hetero,
                         regex_pars = 
                           "^S_act_v_int|^S_act_v_treatconsp|^S_act_v_treathetero|^d_S_fe_av_ncc|^d_S_fe_av_nch|d_S_fe_av_ch|^S_exp_v_int|^S_exp_v_treatconsp|^S_exp_v_treathetero|^d_S_fe_ev_ncc|^d_S_fe_ev_nch|^d_S_fe_ev_ch|^S_pred_v_int|^S_pred_v_treatconsp|^S_pred_v_treathetero|^d_S_fe_pv_ncc|^d_S_fe_pv_nch|^d_S_fe_pv_ch",
                         prob = 0.8, # 80% intervals
                         prob_outer = 0.95, # 95%
                         point_est = "median"
    ) %>% mutate(
      trait = c("Activity", "Exporation", "Predation", "Activity", "Activity", "Exporation", "Exporation", "Predation", "Predation", "Activity", "Exporation", "Predation", "Activity", "Exporation", "Predation", "Activity", "Exporation", "Predation"),
      treat = c("6.No Competition", "6.No Competition", "6.No Competition", "5.Conspecific", "4.Heterospecific", "5.Conspecific", "4.Heterospecific", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "3.NC-Cdiff", "3.NC-Cdiff", "2.NC-Hdiff", "2.NC-Hdiff", "2.NC-Hdiff", "1.C-Hdiff", "1.C-Hdiff", "1.C-Hdiff")
    )
    ),
    aes(y = treat, color = treat)
  ) +
  geom_point(aes(x = m), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = ll, xmax = hh)) +
  geom_linerange(aes(xmin = l, xmax = h), lwd = 2) +
  labs(
    title = bquote("Variance Within Individuals"),
    x = " ", # " \n Posterior Estimate with 80% (thin) 95% (thick) Credible Intervals",
    y = "Treatment" #"Treatment"
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=12))
plot_ACTEXPPRED_feSigmaDiff_nocomp_consp_hetero

#create plot of the repeatability in activity, exloratin and predation assays 
#plot include estimate of the repeatability in activity, exloratin and predation risk taking between treatments


plot_ACTEXPPRED_VarREPDiff_nocomp_consp_hetero <- 
  ggplot(
    (mcmc_intervals_data(mv_act_exp_pred_s_nocomp_consp_hetero_R,
                         regex_pars = 
                           "^R_act_m_treatnocomp|^R_exp_m_treatnocomp|^R_pred_m_treatnocomp|d_R_id_am_ncc|d_R_id_am_nch|d_R_id_am_ch|^R_act_m_treatconsp|^R_exp_m_treatconsp|^R_pred_m_treatconsp|d_R_id_em_ncc|d_R_id_em_nch|d_R_id_em_ch|^R_act_m_treathetero|^R_exp_m_treathetero|^R_pred_m_treathetero|d_R_id_pm_ncc|d_R_id_pm_nch|d_R_id_pm_ch",
                         prob = 0.8, # 80% intervals
                         prob_outer = 0.95, # 95%
                         point_est = "median",
    ) %>% mutate(
      trait = c("Activity", "Activity", "Activity", "Exploration", "Exploration", "Exploration", "Predation", "Predation", "Predation", "Activity", "Exploration", "Predation", "Activity", "Exploration", "Predation", "Activity", "Exploration", "Predation"),
      treat = c("6.No Competition", "5.Conspecific", "4.Heterospecific", "6.No Competition", "5.Conspecific", "4.Heterospecific", "6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "3.NC-Cdiff", "3.NC-Cdiff", "2.NC-Hdiff", "2.NC-Hdiff", "2.NC-Hdiff", "1.C-Hdiff", "1.C-Hdiff", "1.C-Hdiff"),
    )
    ),
    aes(y = treat, color = treat)
  ) +
  geom_point(aes(x = m), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = ll, xmax = hh)) +
  geom_linerange(aes(xmin = l, xmax = h), lwd = 2) +
  labs(
    title = bquote("Behavioural Repeatability "),
    x = "\n Posterior Estimate \n with 80% (thin) 95% (thick) Credible Intervals", # " \n Posterior Estimate with 80% (thin) 95% (thick) Credible Intervals",
    y = " " #"Treatment"
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = 0.001, linetype = "dotted") +
  theme_catplots()
plot_ACTEXPPRED_VarREPDiff_nocomp_consp_hetero

#combine all three plots into one single plot
diff_among_within_rep <- grid.arrange(plot_ACTEXPPRED_VarMeanDiff_nocomp_consp_hetero, plot_ACTEXPPRED_feSigmaDiff_nocomp_consp_hetero, plot_ACTEXPPRED_VarREPDiff_nocomp_consp_hetero,  nrow = 3)
#ggsave("diff_among_within_rep.png", diff_among_within_rep, height = 13, width = 10)

#### fixed effect plots #####

#create plot of overall mean behaviour in activity, exloratin and predation assays 
#plot include estimate of the difference in overall mean behaviour in activity, exloratin and predation risk taking between treatments


plot_ACTEXPPRED_feMeanDiff_nocomp_consp_hetero <- 
  ggplot(
    (mcmc_intervals_data(mv_act_exp_pred_m_nocomp_consp_hetero,
                         regex_pars = 
                           "^M_act_m_int|^M_act_m_treatconsp|^M_act_m_treathetero|^d_M_fe_am_ncc|^d_M_fe_am_nch|d_M_fe_am_ch|^M_exp_m_int|^M_exp_m_treatconsp|^M_exp_m_treathetero|^d_M_fe_em_ncc|^d_M_fe_em_nch|^d_M_fe_em_ch|^M_pred_m_int|^M_pred_m_treatconsp|^M_pred_m_treathetero|^d_M_fe_pm_ncc|^d_M_fe_pm_nch|^d_M_fe_pm_ch",
                         prob = 0.8, # 80% intervals
                         prob_outer = 0.95, # 95%
                         point_est = "median"
    ) %>% mutate(
      trait = c("Activity", "Exporation", "Predation", "Activity", "Activity", "Exporation", "Exporation", "Predation", "Predation", "Activity", "Exporation", "Predation", "Activity", "Exporation", "Predation", "Activity", "Exporation", "Predation"),
      #rep(c("Activity", "Exporation", "Predation"), 6),
      treat = c("6.No Competition", "6.No Competition", "6.No Competition", "5.Conspecific", "4.Heterospecific", "5.Conspecific", "4.Heterospecific", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "3.NC-Cdiff", "3.NC-Cdiff", "2.NC-Hdiff", "2.NC-Hdiff", "2.NC-Hdiff", "1.C-Hdiff", "1.C-Hdiff", "1.C-Hdiff")
      #(rep(c("6.No Competition"), each = 3), rep(c("5.Conspecific", "4.Heterospecific"), times = 3), rep(c("3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"), each = 3))
      #, "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"), each = 3)
    )
    ),
    aes(y = treat, color = treat)
  ) +
  geom_point(aes(x = m), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = ll, xmax = hh)) +
  geom_linerange(aes(xmin = l, xmax = h), lwd = 2) +
  labs(
    title = bquote("Overall Mean Behaviour"),
    x = " ", # " \n Posterior Estimate with 80% (thin) 95% (thick) Credible Intervals",
    y = " " #"Treatment"
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=12))
plot_ACTEXPPRED_feMeanDiff_nocomp_consp_hetero

#create plot of probability of Remaining in the acclimation zone for exloratin and predation assays (hurdle model)
#plot include estimate of the difference in hurdle model for exloratin and predation risk taking between treatments. 

plot_ACTEXPPRED_feHurdleDiff_nocomp_consp_hetero <- 
  ggplot(
    (mcmc_intervals_data(mv_act_exp_pred_m_nocomp_consp_hetero,
                         regex_pars = 
                           "^H_exp_h_int|^H_exp_h_treatconsp|^H_exp_h_treathetero|^d_H_fe_eh_ncc|^d_H_fe_eh_nch|^d_H_fe_eh_ch|^H_pred_h_int|^H_pred_h_treatconsp|^H_pred_h_treathetero|^d_H_fe_ph_ncc|^d_H_fe_ph_nch|^d_H_fe_ph_ch",
                         prob = 0.8, # 80% intervals
                         prob_outer = 0.95, # 95%
                         point_est = "median"
    ) %>% mutate(
      trait = c("Exporation", "Predation", "Exporation", "Exporation", "Predation", "Predation", "Exporation", "Predation", "Exporation", "Predation", "Exporation", "Predation"),
      treat = c("6.No Competition", "6.No Competition", "5.Conspecific", "4.Heterospecific", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "3.NC-Cdiff", "2.NC-Hdiff", "2.NC-Hdiff", "1.C-Hdiff", "1.C-Hdiff")
    )
    ),
    aes(y = treat, color = treat)
  ) +
  geom_point(aes(x = m), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = ll, xmax = hh)) +
  geom_linerange(aes(xmin = l, xmax = h), lwd = 2) +
  labs(
    title = bquote("Probability of Remaining in \n Acclimation Zone (Hurdle model)"),
    x = "\n Posterior Estimate \n with 80% (thin) 95% (thick) Credible Intervals", # " \n Posterior Estimate with 80% (thin) 95% (thick) Credible Intervals",
    y = " " #"Treatment"
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Heterospecific"), bquote(Delta ~ "No Competition-\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=12))
plot_ACTEXPPRED_feHurdleDiff_nocomp_consp_hetero

#combine both plots together
diff_meansFE_mean_hurdle <- grid.arrange(plot_ACTEXPPRED_feMeanDiff_nocomp_consp_hetero, plot_ACTEXPPRED_feHurdleDiff_nocomp_consp_hetero,  nrow = 2)
#ggsave("diff_meansFE_mean_hurdle.png", diff_meansFE_mean_hurdle, height = 13, width = 10)

                                          
