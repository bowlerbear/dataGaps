library(spatstat)
library(raster)
library(tidyverse)
library(tmap)

# Truth ------------------------------------------------------------------------

getTrue <- function(intercept = 6, beta.x = 0, beta.y = 2){
  
bbox <- owin(c(-1,1),c(-1,1))

# Linear predictor of intensity - affected by y
lp <-  function(x,y){ exp(intercept + beta.x * x +  beta.y * y)}

# Point process intensity - points per unit area
lambda <- as.im(lp, bbox)
#plot(lambda)

# Simulate point pattern - total number of points is lambda x area of window
pp <- rpoispp(lambda)
#plot(pp)

# Grid size
n_grids = 20

# Calc points per grid 
grid_raster <- raster(pixellate(pp, dimyx = c(n_grids, n_grids)))
names(grid_raster) <- "state"
grid_df <- as.data.frame(grid_raster, xy=T)
grid_df$site<- 1:nrow(grid_df)
#mean(grid_df$state)

return(grid_df)

}


# Add proxy --------------------------------------------------------------------

#create a degraded y variable

#code taken from:
#https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
complement <- function(y, rho, x) {
  if (missing(x)) x <- runif(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

addProxy <- function(grid_df, rho = 0.7){

grid_df$proxy_y <- complement(y = grid_df$y, rho = rho)
#cor(grid_df$y, grid_df$proxy_y)

#cut into 20 values to mirror the current y
grid_df$proxy_y <- as.numeric(cut(grid_df$proxy_y, 20))

return(grid_df)

}

# Data -------------------------------------------------------------------------

getMAR <- function(grid_df, intercept = -1, beta.x = 0, beta.y = 4, size = 100){ 

# Probability of grid cell being sampled - also affected by y
grid_df$prob_sampled <- boot::inv.logit(intercept + beta.x * grid_df$x + beta.y * grid_df$y) 
#qplot(y, x,data=grid_df, colour=prob_sampled)
#means we are more likely to sample higher values

#take a sample
sampled_sites <- sample(grid_df$site, size = size, replace = FALSE, prob = grid_df$prob_sampled)
grid_df$sampled <- ifelse(grid_df$site %in% sampled_sites, 1, 0)

# Insert missing values into unsampled cells
grid_df_mar <- grid_df %>% mutate(count = ifelse(sampled==0, NA, state))

return(grid_df_mar)

}

# Jags model code --------------------------------------------------------------

cat("model{
  for(i in 1:N){ 

#data distribution - count has NAs
count[i] ~ dpois(lambda[i])

#linear predictors - impute
log(lambda[i]) <- alpha.lam + beta.psi * y[i]

}

#get overall mean
meanAbund <- mean(lambda)
 
# priors 
alpha.lam ~ dnorm(0, 0.01) 
beta.psi ~ dnorm(0, 0.01) 

}",file="imputeModel.txt")


# Solutions -------------------------------------------------------------------- 

getModels <- function(grid_df, grid_df_mar){
  
grid_df_mar_sampled <- grid_df_mar %>% filter(sampled==1) 

# true state -------------------------------------------------------------------

state_df = data.frame(type = "true",
                      mean = mean(grid_df$state),
                      lower = NA,
                      upper = NA)

# naive ------------------------------------------------------------------------

mean(grid_df$state)
mean(grid_df_mar_sampled$count)

naive_glm1 <- glm(count ~ 1, data = grid_df_mar_sampled, family = poisson)
exp(confint(naive_glm1))

naive_df = data.frame(type = "naive",
                      mean = exp(coef(naive_glm1)),
                      lower = exp(confint(naive_glm1))[1],
                      upper = exp(confint(naive_glm1))[2])

# subsampling ------------------------------------------------------------------

#subsampling the data according to y

#ideally we want at least 1 data points for each y value - not fully possible

subsample_fun <- function(grid_df_mar_sampled){
  
  grid_df_mar_sampled %>%
    group_by(y) %>% #the axis that affects biodiversity
    #sample_n(size = ifelse(length(count)>=y_25, y_25, length(count))) %>%
    sample_n(size = 1) %>%
    ungroup() %>%
    summarise(mean_count = mean(count)) %>%
    pull(mean_count)
  
}

#hist(grid_df_mar_subsampled$y)
#mean(grid_df_mar_subsampled$count)

#apply function
grid_df_mar_subsampled <- replicate(100, subsample_fun(grid_df_mar_sampled)) 

subsampled_df = data.frame(type = "subsampled",
                           mean = mean(grid_df_mar_subsampled),
                           lower = quantile(grid_df_mar_subsampled, 0.025),
                           upper = quantile(grid_df_mar_subsampled, 0.975))

# weighting ---------------------------------------------------------------------

#model weights

#hist(1/grid_df_mar_sampled$prob_sampled)
weighted_glm1 <- glm(count ~ 1, data = grid_df_mar_sampled, family = poisson,
                     weights = 1/prob_sampled)

#use sandwich estimators for SE
exp(coefci(weighted_glm1, vcov = sandwich))

weighted_df = data.frame(type = "weighted",
                         mean = exp(coef(weighted_glm1)),
                         lower = exp(coefci(weighted_glm1, vcov = sandwich))[1],
                         upper = exp(coefci(weighted_glm1, vcov = sandwich))[2])

#quasi-random

weighted_design <- svydesign(ids=~1,
                             data = grid_df_mar_sampled,
                             probs = ~ prob_sampled)

svyglm1 <- svyglm(count ~ 1, design = weighted_design, family = poisson)
exp(confint(svyglm1))
#this gives the same as above

svy_df = data.frame(type = "svyglm",
                   mean = exp(coef(svyglm1)),
                   lower = exp(confint(svyglm1))[1],
                   upper = exp(confint(svyglm1))[2])

#post stratify

basic_design <- svydesign(ids = ~1,
                          data =  grid_df_mar_sampled)

ps.weights <- data.frame(y = names(table(grid_df_mar$y)),
                         Freq = as.numeric(table(grid_df_mar$y)))

ps_design <- postStratify(design = basic_design,
                          strata = ~ y,
                          population = ps.weights,
                          partial = TRUE)

(ps_samp_mean <- svymean(design = ps_design, x=~count))
confint(object = ps_samp_mean, level = 0.95)

ps_df = data.frame(type = "poststratified",
                   mean = ps_samp_mean[1],
                   lower = confint(object = ps_samp_mean, level = 0.95)[1],
                   upper = confint(object = ps_samp_mean, level = 0.95)[2])

# imputation -------------------------------------------------------------------

bugs.data <- list(count = grid_df_mar$count,
                  N = nrow(grid_df_mar),
                  y = grid_df_mar$y)

jags.model <- jags(bugs.data, inits=NULL, "meanAbund", 
                   model.file = "imputeModel.txt", 
                   n.chains = 3, n.iter = 500, verbose=FALSE)

imputed_df = data.frame(type = "imputed",
                        mean = jags.model$summary["meanAbund",1],
                        lower = jags.model$summary["meanAbund",3],
                        upper = jags.model$summary["meanAbund",7])

# combine all ------------------------------------------------------------------  

results_df <- rbind(state_df,
                    naive_df,
                    subsampled_df,
                    weighted_df,
                    svy_df,
                    ps_df,
                    imputed_df)

return(results_df)

}





