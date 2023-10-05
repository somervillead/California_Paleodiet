### MixSIAR Model for southern California stable isotope data

### Andrew D. Somerville
### Using r package, MixSIar 

### MixSIAR developed by BC Stock, and explained in:
### Stock, B. C., Jackson, A. L., Ward, E. J., Parnell, A. C., Phillips, D. L., & Semmens, B. X. (2018). Analyzing mixing systems using a new generation of Bayesian tracer mixing models. PeerJ, 6, e5096.
### Stock, B. C. and B. X. Semmens (2016). MixSIAR GUI User Manual.Version 3.1. https://github.com/brianstock/MixSIAR/. doi:10.5281/zenodo.47719.


### This is the full mixing model that includes inland, coastal, and island sites

library(tidyverse)
library(MixSIAR)


# Loading mix data

mix <- load_mix_data(filename="CONSUMER.csv",
                     iso_names=c("d13C","d15N"),
                     factors=c("Prov", "Period"),
                     fac_random=c(TRUE, TRUE),
                     fac_nested=c(FALSE, TRUE),
                     cont_effects=NULL)

# Loading source data

source <- load_source_data(filename="SOURCE_20201012.csv",
                           source_factors=NULL,
                           conc_dep=TRUE,
                           data_type="means",
                           mix)

# Loading discrimination/TDF data

discr <- load_discr_data(filename="TEF.csv", mix)

# Isospace plot

plot_data(filename="isospace_plot",
          plot_save_pdf=TRUE,
          plot_save_png=FALSE,
          mix,source,discr)


################################################################################
# # INFORMATIVE prior (construct alpha from geographic assumptions [island, coastal, inland])
################################################################################

# Resource importance scored as high = 3, med = 2, low = 1
# In the individual models for island, coastal, and inland sites, informed priors were used for each setting. For this general model, however, I use uninformed priors, ie. 1,1,1,1.


# Uninformed Prior
Cali.alpha <- c(1,1,1,1)

# Generate alpha hyperparameters scaling sum(alpha)=n.sources
Cali.alpha <- Cali.alpha*length(Cali.alpha)/sum(Cali.alpha)

# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)
Cali.alpha[which(Cali.alpha==0)] <- 0.01

# Plot informative prior
plot_prior(alpha.prior=Cali.alpha,
           source=source,
           plot_save_pdf=FALSE,
           plot_save_png=FALSE,
           filename="prior_plot_Cali_inf")


# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- FALSE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("test" first, then "short", then "normal")
#jags.1 <- run_model(run="normal", mix, source, discr, model_filename, alpha.prior=1)
#jags.1 <- run_model(run="short", mix, source, discr, model_filename, alpha.prior=1)
jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior=1)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.1, mix, source)


###################################################

## OUTPUT PLOTTING ##


library(R2jags)
library(tidyverse)


attach.jags(jags.1)
post.Coastal <- data.frame(Location = "Coastal", Marine.High = p.fac1[,1,1], Marine.Low = p.fac1[,1,2], Plants = p.fac1[,1,3], T.Mammal = p.fac1[,1,4])
post.Inland <- data.frame(Location = "Inland", Marine.High = p.fac1[,2,1], Marine.Low = p.fac1[,2,2], Plants = p.fac1[,2,3], T.Mammal = p.fac1[,2,4])
post.Island <- data.frame(Location = "Island", Marine.High = p.fac1[,3,1], Marine.Low = p.fac1[,3,2], Plants = p.fac1[,3,3], T.Mammal = p.fac1[,3,4])


Coastal <- post.Coastal %>% gather(source,value, 2:5)
Inland <-post.Inland %>% gather(source,value, 2:5)
Island <- post.Island %>% gather(source,value, 2:5)


all <- rbind(Coastal, Inland, Island)


########## OVERALL VIOLIN PLOTS ###########

ggplot(aes(x = value, y=source, fill=source),data=all)+ 
  geom_violin(trim = FALSE,
        alpha = 0.8)+
  geom_boxplot(outlier.colour=NA, alpha=0.3, color="black", width=.1)+
    theme_bw() +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position="none", legend.text=element_text(size=16), legend.key.size = unit(1, "cm"))+
  scale_fill_manual(values=c("#225ea8", "#41b6c4", "#7fcdbb", "#ffffd9"))+ 
  facet_wrap(~Location, ncol=2, scales="free") +
  coord_flip()+
  xlab("Diet Proportions") + 
  ggsave(filename = "Cali_Violin_Plots_CONSUMER_noSNI__LOCATION_violin.tiff", width=8, height=7)


