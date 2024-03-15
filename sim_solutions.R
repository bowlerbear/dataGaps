# libraries --------------------------------------------------------------------

require(survey)
library(rjags)
library(jagsUI)
library(sandwich)
library(lmtest)

# function file ----------------------------------------------------------------

source("sim_functions.R")

# get truths -------------------------------------------------------------------

#for each simulation run, generate a community

truth_df <- lapply(1:100, function(i){ 
                getTrue(intercept = 6, beta.x = 0, beta.y = 2) %>%
                  add_column(nu = i)
}) %>% bind_rows()

#change to beta.x = 1 for the SI plot with an additional covariate

# parameter set up -------------------------------------------------------------
        
# to be applied for each truth          
params <- expand.grid(sample_prop = seq(0.05,0.95, 0.1),
                      proxy_corr = seq(0.05, 0.95, 0.1),
                      sim_nu = sort(unique(truth_df$nu)))
                  
# get data ---------------------------------------------------------------------

all_results <- lapply(1:nrow(params), function(i){

#params for this sim
sim_sample_prop = params$sample_prop[i]
sim_proxy_corr = params$proxy_corr[i]  
sim_nu = params$sim_nu[i] 
  
#get truth version for this run
grid_df <- truth_df %>%
              filter(nu == sim_nu)
  
# sampling and abundance affected by y
grid_df_mar <- getMAR(grid_df, intercept = 0.2, beta.x = 0, beta.y = 2, 
                      size = round(nrow(grid_df) * sim_sample_prop))

#proxy correlation if we are using an imperfect variable
if(sim_proxy_corr<1){
  grid_df_mar <- addProxy(grid_df_mar, rho = sim_proxy_corr)
  grid_df_mar$y <- grid_df_mar$proxy_y
}

#get post hoc inclusion probabilities
sampled_glm <- glm(sampled ~ y, data = grid_df_mar, family="binomial")
grid_df_mar$prob_sampled <- predict(sampled_glm, type="response")

# get attempted solutions
results_df <- getModels(grid_df, grid_df_mar)

#add on param values
results_df %>%
  add_column(proxy_corr = sim_proxy_corr,
             sample_prop = sim_sample_prop,
             sim_nu = sim_nu)

}) %>% bind_rows()

# summarising ------------------------------------------------------------------

all_results_summary <- all_results %>%
  group_by(sim_nu, proxy_corr, sample_prop) %>%
  mutate(truth = mean[type=="true"]) %>%
  group_by(type, proxy_corr, sample_prop) %>%
  mutate(ci_width = upper - lower,
         bias = mean - truth) %>%
  summarise(mean = median(bias),
            ci_mean = median(ci_width)) %>%
  ungroup() %>%
  filter(type!="true") %>%
  filter(type!="svyglm") %>% # same as weighted
  mutate(type = factor(type, levels=c("naive","subsampled","weighted","poststratified","imputed")))

# all facets - mean  --------------------------------------------------------

all_results_summary %>%
  ggplot(aes(x = sample_prop, y = mean,colour = type)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_wrap(~ proxy_corr) +
  ylab("Predicted mean abundance") +
  xlab("Sampled fraction")

# all facets - uncertainty -------------------------------------------------

all_results_summary %>%
  ggplot(aes(x = sample_prop, y = ci_mean, colour = type)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_wrap(~ proxy_corr) +
  ylab("CI width of predicted mean abundance") +
  xlab("Sampled fraction")

# simple plots -----------------------------------------------------------------------

all_results_summary$type <- as.character(all_results_summary$type)
all_results_summary$type[all_results_summary$type=="weighted"] <- "weighted glm"
all_results_summary$type <- factor(all_results_summary$type,
                                   levels=c("naive","subsampled","weighted glm",
                                            "poststratified","imputed"))
g1 <- all_results_summary %>%
  filter(proxy_corr>0.7 & proxy_corr<0.8) %>%
  ggplot(aes(x = sample_prop, y = mean, colour = type)) +
  geom_point() +
  geom_line() +
  scale_color_brewer("Model",palette= "Set2")+
  theme_classic() +
  ylab("Bias (deviation from truth)") +
  xlab("Sampled fraction")

g2 <- all_results_summary %>%
  filter(sample_prop>0.3 & sample_prop<0.4) %>%
  ggplot(aes(x = proxy_corr, y = mean, colour = type)) +
  geom_point() +
  geom_line() +
  scale_color_brewer("Model",palette= "Set2")+
  theme_classic() +
  ylab("Bias (deviation from truth)") +
  xlab("Covariate knowledge (correlation with truth)")

meanPlot <- cowplot::plot_grid(g1,g2, labels=c("A","C"), nrow=2)

# uncertainty -----------------------------------------------------------------------

g1 <- all_results_summary %>%
  filter(proxy_corr>0.7 & proxy_corr<0.8) %>%
  ggplot(aes(x = sample_prop, y = ci_mean, colour = type)) +
  geom_point() +
  geom_line() +
  scale_color_brewer("Model",palette= "Set2")+
  theme_classic() +
  ylab("Uncertainty (CI width)") +
  xlab("Sampled fraction")

g2 <- all_results_summary %>%
  filter(sample_prop>0.3 & sample_prop<0.4) %>%
  ggplot(aes(x = proxy_corr, y = ci_mean, colour = type)) +
  geom_point() +
  geom_line() +
  scale_color_brewer("Model",palette= "Set2")+
  theme_classic() +
  ylab("Uncertainty (CI width)") +
  xlab("Covariate knowledge (correlation with truth)")

ciPlot <- cowplot::plot_grid(g1,g2, labels=c("B","D"), nrow=2)

# combining -------------------------------------------------------

cowplot::plot_grid(meanPlot,
                   ciPlot,
                   ncol = 2)

ggsave("plots/Fig.5.png",width = 11, height = 6)
# end --------------------------------------------------------------------------