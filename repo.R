###############################################################################################
###################################### Brian Forester #########################################
################# Globalization and the Demand for Joint Military Exercises ###################
###############################################################################################

## READ ME: This script reproduces the results of "Globalization and the Demand for Joint Military Exercises." 
## This is currently an unpublished manuscript. Please contact the author for manuscript access.
## Files needed for reproduction: ctryYrData.rda, functions.R, repo.R

# Clear the environment
rm(list=ls())

# Set random seed
seed <- 99
set.seed(seed)

# Set working directory. 
dir='INSERT YOUR WORKING DIRECTORY HERE'

# Load standard dependencies and helpers from the functions.R script
source(paste0(dir, 'functions.R'))

# Additional dependencies required
loadPkg(c('sandwich', 'lmtest',             # VCOV adjustments
          'MASS',                           # Negative binomial models
          'lme4'))                          # Random effects      

# Load the country year data, 1978-2016
load(file=paste0(dir,'ctryYrData.rda'))


#################### Figure 1 -- JME participation, by regime ###################


# Organize data for plotting
plotData <- ctryYrData %>% 
        filter(ccode != 2, ccode != 200, ccode != 220,
               ccode != 365, ccode != 710, 
               ccode != 255, ccode != 740) %>%
        group_by(year) %>%
        summarise(dem_jme=sum(jme_d[democracy_vdem==1], na.rm=T),
                  dem_nojme=dem_jme + sum(democracy_vdem[jme_d==0], na.rm=T),
                  autoc_jme=dem_nojme + sum(jme_d[autoc==1], na.rm=T),
                  autoc_nojme=autoc_jme + sum(autoc[jme_d==0], na.rm=T))

# JME ribbon plot
jmeribbon <- ggplot(plotData, aes(x=year)) +
        geom_ribbon(aes(ymin=dem_nojme-dem_jme, ymax=dem_nojme), fill='gray30', alpha=.4) +
        geom_ribbon(aes(ymin=0, ymax=dem_nojme-dem_jme), fill='gray80', alpha=.4) +
        geom_ribbon(aes(ymin=dem_nojme, ymax=autoc_jme), fill='gray30', alpha=.4) +
        geom_ribbon(aes(ymin=autoc_jme, ymax=autoc_nojme), fill='gray80', alpha=.4) + 
        geom_line(aes(y=dem_nojme), lwd=.75) +
        scale_x_continuous(breaks=seq(1980, 2015, 5)) +
        scale_y_continuous(breaks=seq(0, 160, 40)) +
        geom_text(x=1983, y=120, label="Autocracies", size=3, fontface='bold')+
        geom_text(x=1983, y=10, label="Democracies", size=3, fontface='bold')+
        theme_classic()+
        labs(x="Year", y="Number of Minor Powers") +
        theme(axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"))

# Export to working directory
ggsave(jmeribbon, file=paste0(dir, 'jmeribbon.pdf'), height=4, width=7)


############### Figure 2 -- KOF over time ##############


## Look at state KOF over time, by autocracies and democracies
plotData <- ctryYrData %>% 
        group_by(year) %>%  
        summarise(kof_a = mean(KOF_econ[democracy_vdem==0], na.rm=T),
                  kof_na = mean(KOF_econ[democracy_vdem==1], na.rm=T),
                  kof = mean(KOF_econ, na.rm=T))

# Plot KOF over time, by autoc and non-autoc
kofovertime <- ggplot(plotData, aes(x=year)) +
        geom_jitter(aes(y=kof_a), size=1, alpha=.5, color="darkslategray") +
        geom_jitter(aes(y=kof_na), size=1, alpha=.5, color="darkslategray") +
        geom_smooth(aes(y=kof), size=.5, linetype=2, color="black") +
        geom_text(x=1990, y=55, label="Democracies", size=3, fontface='bold') +
        geom_text(x=2010, y=45, label="Autocracies", size=3, fontface='bold') +
        scale_x_continuous(breaks=seq(1980, 2015, 5)) +
        theme_classic()+
        labs(x="Year", y="Mean KOF Economic Index") +
        theme(axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"))
kofovertime

# Export to working directory
ggsave(kofovertime, file=paste0(dir, 'kofovertime.pdf'), height=4, width=7)


############### Table 1 -- Summary Statistics ##############


# Subset to minor powers for all further analysis
modelData <- ctryYrData %>% filter(mp==0)

# Export summary statistics table
stargazer(as.data.frame(modelData),
          title="Summary Statistics",
          covariate.labels = c("Joint Military Exercise (binary)", 
                               "KOF Economic Index", "Democracy", 
                               "Rivalries (count)", "Alliances",
                               "Mil. Expenditure (log)",
                               "GDPPC (log)"),
          omit=c(1:4,6:24, 26, 28:29, 33:34),
          omit.summary.stat = c("median", "p25", "p75"))


############### Table 2 -- Main Models ##############


# Model 1 -- Baseline
m1 <- glm(jme_d ~ KOF_econ_lag  +  
                  time + time2 + time3,
          data=modelData,family=binomial('logit'))

# Model 2 -- Baseline with democ
m2 <- glm(jme_d ~ KOF_econ_lag +  
                  democracy_vdem_lag + 
                  time + time2 + time3,
          data=modelData,family=binomial('logit'))

# Model 3 -- JME full controls
m3 <- glm(jme_d ~ KOF_econ_lag  +  
                  democracy_vdem_lag + log_gdppc_lag + log_milex_lag + 
                  alliances_lag + rivalries_lag +
                  time + time2 + time3,
          data=modelData,family=binomial('logit'))

# Model 4 -- JME interaction 
m4 <- glm(jme_d ~ KOF_econ_lag  *  
                  democracy_vdem_lag + log_gdppc_lag + log_milex_lag + 
                  alliances_lag + rivalries_lag + 
                  time + time2 + time3,
          data=modelData,family=binomial('logit'))

# Compute clustered SEs
m1SE <- sqrt(diag(vcovCL(m1, cluster=~wto)))
m2SE <- sqrt(diag(vcovCL(m2, cluster=~wto)))
m3SE <- sqrt(diag(vcovCL(m3, cluster=~wto)))
m4SE <- sqrt(diag(vcovCL(m4, cluster=~wto)))

# Export to latex table
stargazer(m1, m2, m3, m4, 
          title = "Minor Power JME Participation, 1978-2016",
          covariate.labels = c("KOF Economic Index", "Democracy", "GDPPC (log)",
                               "Military Expenditures (log)", "Alliances",
                               "Rivalries (count)",  
                               "Time", "Time$^{2}$", "Time$^{3}$",
                               "KOF x Democracy"),
          column.separate = c(1,1,1,1),
          column.labels = c("Base", "Democracy", "Controls", "Interaction"),
          dep.var.caption="",
          dep.var.labels.include = FALSE,
          add.lines=list(c('Countries',161,161,161,161)),
          model.numbers = TRUE,
          se=list(m1SE, m2SE, m3SE, m4SE),
          p.auto=TRUE,
          digits=2,
          star.cutoffs = .05,
          model.names=FALSE,
          omit.stat = c("f", "ser", 'bic'),
          notes=c("$^{*}$p$<$0.05; SEs clustered on 161 countries; all predictors lagged one year."),
          notes.append=FALSE,
          notes.align="l",
          table.placement = "H")


##################### Model viz -- pred. prob. plots ##################


# Vector of KOF quantities of interest
kofVec <- seq(min(modelData$KOF_econ_lag, na.rm=T), 
              max(modelData$KOF_econ_lag, na.rm=T), 
              length.out = 100)

# Function to generate logit predicted probability plots
logitPredPlot <- function(models, data, vec, sims, interaction){
        
        # No interaction models
        if(interaction==FALSE){
                
                # Generate df of pred probs for each model in the list
                plotData <- lapply(seq_along(models), function(i){
                        
                        # draw betas, incorporate robust standard errors
                        betas <- mvrnorm(n=sims, 
                                         mu=coef(models[[i]]), 
                                         Sigma=vcovCL(models[[i]], cluster=~wto))
                        
                        # generate pred. probs using betas
                        p <- plogis(as.matrix(data) %*% t(betas)) 
                        
                        # organize into a df, categorizing by model name.
                        data.frame(vec = vec,
                                   pe=apply(p,1,mean),
                                   lwr=apply(p,1,quantile, prob=.025),
                                   upr=apply(p,1,quantile, prob=.975),
                                   model=names(models)[[i]])
                })
                
                # Extract plotdata into a single df
                plotData <- do.call(rbind, plotData)
                
                # Plot
                plot <- ggplot(plotData, aes(x=vec)) + 
                        geom_ribbon(aes(ymin = lwr, ymax = upr), fill='gray80', alpha = .4) +
                        geom_line(aes(y=pe)) +
                        geom_rug(data=modelData, aes(x=KOF_econ_lag), sides='b', 
                                 alpha=0.05) +
                        labs(x="KOF Economic Index", y="Pr(JME)")+
                        scale_x_continuous(breaks=seq(0, 100, 10)) +
                        theme(plot.title=element_text(hjust=0.5), 
                              axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"),
                              legend.title=element_blank(),
                              legend.position="bottom",
                              panel.grid.minor.x = element_blank(),
                              panel.grid.minor.y = element_blank())
                
                # Export to working directory
                ggsave(plot, file=paste0(dir, length(models), 'nointeraction.pdf'), width=7, height=5)
                plot
        }
        
        # Interaction models
        else if(interaction==TRUE){
                
                
                # Generate df of pred probs for each model in the list
                plotData <- lapply(seq_along(models), function(i){
                        
                        # draw betas
                        betas <- mvrnorm(n=sims, 
                                         mu=coef(models[[i]]), 
                                         Sigma=vcovCL(models[[i]], cluster=~wto))
                        
                        # for each dataset 
                        tmp <- lapply(seq_along(data), function(j){
                                
                                # generate pred. probs using betas
                                p <-  plogis(as.matrix(data[[j]]) %*% t(betas)) 
                                
                                # baseline pred. probs against the minimum, to reflect change in prob.
                                data.frame(vec = vec,
                                           pe=apply(p,1,mean),
                                           lwr=apply(p,1,quantile, prob=.025),
                                           upr=apply(p,1,quantile, prob=.975),
                                           model=names(models)[[i]],
                                           cat=names(data)[[j]])
                        })
                        do.call(rbind, tmp)
                })
                
                # Extract plotdata into a single df
                plotData <- do.call(rbind, plotData)
                
                # Plot
                plot <- ggplot(plotData, aes(x=vec)) + 
                        geom_ribbon(aes(ymin = lwr, ymax = upr, fill=cat), alpha = .4) +
                        geom_line(aes(y=pe, linetype=cat)) +
                        geom_rug(data=modelData, aes(x=KOF_econ_lag), sides='b', 
                                 alpha=0.05) +
                        labs(x="KOF Economic Index", y="Pr(JME)")+
                        scale_x_continuous(breaks=seq(0, 100, 10))+
                        scale_fill_manual(values=c('gray80', 'gray80')) +
                        theme(plot.title=element_text(hjust=0.5), 
                              axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"),
                              legend.title=element_blank(),
                              legend.position="bottom",
                              panel.grid.minor.x = element_blank(),
                              panel.grid.minor.y = element_blank())
                
                # Export to working directory
                ggsave(plot, file=paste0(dir, 'interaction.pdf'), width=7, height=5)
                plot
        }
}


############### Figure 3 -- model 3 #############################


# Sim data for models without the interaction
genData <- with(modelData, data.frame(intercept=1,
                                      KOF_econ_lag=kofVec,
                                      democracy_vdem_lag=0,
                                      log_gdppc_lag=mean(log_gdppc_lag, na.rm=T),
                                      log_milex_lag=mean(log_milex_lag, na.rm=T),
                                      alliances_lag=1,
                                      rivalries_lag=1,
                                      time=mean(time),
                                      time2=mean(time)^2,
                                      time3=mean(time)^3))

# Package model 3 into a list for function input
noInteractionMods <- list(m3)
names(noInteractionMods) <- c('JME')

# Run the function on model 3
logitPredPlot(noInteractionMods, genData, kofVec, sims=1000, interaction=F)


############### Figure 4 -- model 4 (interaction) ##############


# Sim data with interaction -- democracy
demData <- with(modelData, data.frame(intercept=1,
                                      KOF_econ_lag=kofVec,
                                      democracy_vdem_lag=1,
                                      log_gdppc_lag=mean(log_gdppc_lag, na.rm=T),
                                      log_milex_lag=mean(log_milex_lag, na.rm=T),
                                      alliances_lag=1,
                                      rivalries_lag=1,
                                      time=mean(time),
                                      time2=mean(time)^2,
                                      time3=mean(time)^3,
                                      interaction=kofVec))

# Sim data with interaction -- autocracy
nonDemData <- with(modelData, data.frame(intercept=1,
                                         KOF_econ_lag=kofVec,
                                         democracy_vdem_lag=0,
                                         log_gdppc_lag=mean(log_gdppc_lag, na.rm=T),
                                         log_milex_lag=mean(log_milex_lag, na.rm=T),
                                         alliances_lag=1,
                                         rivalries_lag=1,
                                         time=mean(time),
                                         time2=mean(time)^2,
                                         time3=mean(time)^3,
                                         interaction=0))

# Package model 4 as a list for function input
interactionMods <- list(m4)
names(interactionMods) <- c('JME')

# Package the interaction data into a list
InteractionData <- list(demData, nonDemData)
names(InteractionData) <- c('Democracy', 'Autocracy')

# Run the function on model 4
logitPredPlot(interactionMods, InteractionData, kofVec, sims=1000, interaction=T)


############### Supplemental Information -- Robustness Checks ##############


#### Table 3 -- Random Effects  Models ####

# Model 1 -- Baseline
m1RE <- glmer(jme_d ~ scale(KOF_econ_lag)  +  
                      scale(time) + scale(time2) + scale(time3) + (1 | wto),
              control=glmerControl('bobyqa', optCtrl=list(maxfun=2e5)),
              data=modelData,family=binomial('logit'))

# Model 2 -- with Democ
m2RE <- glmer(jme_d ~ scale(KOF_econ_lag)  +  
                      democracy_vdem_lag + 
                      scale(time) + scale(time2) + scale(time3) + (1 | wto),
              control=glmerControl('bobyqa', optCtrl=list(maxfun=2e5)),
              data=modelData,family=binomial('logit'))

# Model 3 -- Full controls
m3RE <- glmer(jme_d ~ scale(KOF_econ_lag) +
                      democracy_vdem_lag + scale(log_gdppc_lag) + scale(log_milex_lag) + 
                      alliances_lag + rivalries_lag +  
                      scale(time) + scale(time2) + scale(time3) + (1 | wto),
              control=glmerControl('bobyqa', optCtrl=list(maxfun=2e5)),
              data=modelData,family=binomial('logit'))

# Model 4 --  interaction
m4RE <- glmer(jme_d ~ scale(KOF_econ_lag) * democracy_vdem_lag +  
                      democracy_vdem_lag + scale(log_gdppc_lag) + scale(log_milex_lag) + 
                      alliances_lag + rivalries_lag +  
                      scale(time) + scale(time2) + scale(time3) + (1 | wto),
              control=glmerControl('bobyqa', optCtrl=list(maxfun=2e5)),
              data=modelData,family=binomial('logit'))

# Table 3 latex
stargazer(m1RE, m2RE, m3RE, m4RE,
          title = "Model Specifications with Country Random Effects",
          covariate.labels = c("KOF Economic Index", "Democracy", "GDPPC (log)",
                               "Military Expenditures (log)", "Alliances",
                               "Rivalries (count)",  
                               "Time", "Time$^{2}$", "Time$^{3}$",
                               "KOF x Democracy"),
          column.separate = c(1,1,1,1),
          column.labels = c("Base", "Democracy", "Controls", "Interaction"),
          dep.var.caption="",
          dep.var.labels.include = FALSE,
          add.lines=list(c('Countries',161,161,161,161),
                         c('sigma (country)', round(m1RE@theta,3), round(m2RE@theta,3), 
                           round(m3RE@theta,3), round(m4RE@theta,3))), #RE variance
          model.numbers = FALSE,
          p.auto=TRUE,
          digits=3,
          star.cutoffs = .05,
          model.names=FALSE,
          omit.stat = c("f", "ser", 'bic'),
          notes=c("$^{*}$p$<$0.05; continuous predictors scaled; all predictors lagged one year."),
          notes.append=FALSE,
          notes.align="l",
          table.placement = "H")


#### Table 4 -- Count Models -- Negative Binomial ####

# Negative binomial -- M1 base
negBin1 <- glm.nb(data = modelData, count ~  KOF_econ_lag  + as.factor(year),
                  link=log)

# Negative binomial -- M2 democ
negBin2 <- glm.nb(data = modelData, count ~  KOF_econ_lag  + 
                          democracy_vdem_lag  + 
                          as.factor(year),
                  link=log)

# Negative binomial -- M3 full
negBin3 <- glm.nb(data = modelData, count ~  KOF_econ_lag  + 
                          democracy_vdem_lag  + log_gdppc_lag + 
                          log_milex_lag + 
                          alliances_lag + rivalries_lag +
                          as.factor(year),
                  link=log)

# Negative binomial -- M4 interaction
negBin4 <- glm.nb(data = modelData, count ~  KOF_econ_lag  *  
                          democracy_vdem_lag + log_gdppc_lag + 
                          log_milex_lag + 
                          alliances_lag + rivalries_lag +  as.factor(year),
                  link=log)

# Extract coefficients, excluding the time FEs
negBin1Coef <- coef(negBin1)[1:2]
negBin2Coef <- coef(negBin2)[1:3]
negBin3Coef <- coef(negBin3)[1:7]
negBin4Coef <- coef(negBin4)[c(1:7, length(coef(negBin4)))] # interaction term

# Compute clustered standard errors
negBin1SE <- sqrt(diag(vcovCL(negBin1, cluster=~wto)))[1:2]
negBin2SE <- sqrt(diag(vcovCL(negBin2, cluster=~wto)))[1:3]
negBin3SE <- sqrt(diag(vcovCL(negBin3, cluster=~wto)))[1:7]
negBin4SE <- sqrt(diag(vcovCL(negBin4, cluster=~wto)))[c(1:7, length(coef(negBin4)))] # interaction term

# Table 4 latex
stargazer(negBin1, negBin2, negBin3, negBin4,
          title = "Count Models of JME Participation",
          covariate.labels = c("KOF Economic Index", "Democracy", "GDPPC (log)", 
                               "Military Expenditures (log)",  "Alliances",
                               "Rivalries (count)",  
                               "KOF x Democracy"),
          column.separate = c(1, 1, 1, 1),
          column.labels = c("Base", "Democracy", "Controls", "Interaction"),
          dep.var.caption="",
          dep.var.labels.include = FALSE,
          add.lines=list(c('Countries',161,161,161,161)),
          model.numbers = FALSE,
          coef=list(negBin1Coef, negBin2Coef, negBin3Coef, negBin4Coef),
          se=list(negBin1SE, negBin2SE, negBin3SE, negBin4SE),
          omit=c('as\\.factor\\(year\\)', 'time'), #regex that omits fixed effects
          p.auto=TRUE,
          digits=3,
          star.cutoffs = .05,
          model.names=FALSE,
          omit.stat = c("f", "ser", 'bic'),
          notes=c("$^{*}$p$<$0.05; Year FEs not depicted; SEs clustered on 161 countries; all predictors lagged one year."),
          notes.append=FALSE,
          notes.align="l",
          table.placement = "H")


#### Table 5 -- Other specifications to assess robustness ####

# Continuous measure of democracy
contDem <- glm(jme_d ~ KOF_econ_lag  +  
                       v2x_polyarchy_lag + log_gdppc_lag + log_milex_lag + 
                       alliances_lag + rivalries_lag +  
                       time + time2 + time3,
               data=modelData,family=binomial('logit'))

# Polity measure of democ
polityDem <- glm(jme_d ~ KOF_econ_lag  +  
                         democracy_polity_lag + log_gdppc_lag + log_milex_lag + 
                         alliances_lag + rivalries_lag +  
                         time + time2 + time3,
                 data=modelData,family=binomial('logit'))

# JME Lag included
jmeLag <- glm(jme_d ~ KOF_econ_lag  + jme_d_lag +  
                      democracy_vdem_lag + log_gdppc_lag + log_milex_lag + 
                      alliances_lag + rivalries_lag +  
                      time + time2 + time3,
              data=modelData,family=binomial('logit'))

# CINC in lieu of milex
cinc <- glm(jme_d ~ KOF_econ_lag  +  
                    democracy_vdem_lag + log_gdppc_lag + cinc_lag + 
                    alliances_lag + rivalries_lag +  
                    time + time2 + time3,
            data=modelData,family=binomial('logit'))

# GDP in lieu of GDPPC
gdp <- glm(jme_d ~ KOF_econ_lag  +  
                   democracy_vdem_lag + log_gdp_lag + log_milex_lag + 
                   alliances_lag + rivalries_lag +  
                   time + time2 + time3,
           data=modelData,family=binomial('logit'))

# Extract coefficients
contDemCoef <- coef(contDem)
polityDemCoef <- coef(polityDem)
jmeLagCoef <- coef(jmeLag)
gdpCoef <- coef(gdp)
cincCoef <- coef(cinc)

# Compute clustered standard errors
contDemSE <- sqrt(diag(vcovCL(contDem, cluster=~wto)))
polityDemSE <- sqrt(diag(vcovCL(polityDem, cluster=~wto)))
jmeLagSE <- sqrt(diag(vcovCL(jmeLag, cluster=~wto)))
gdpSE <- sqrt(diag(vcovCL(gdp, cluster=~wto)))
cincSE <- sqrt(diag(vcovCL(cinc, cluster=~wto)))


# Table 5 latex
stargazer(contDem, polityDem, jmeLag, gdp, cinc,
          title = "Alternative Model Specifications",
          covariate.labels = c("KOF Economic Index",
                               "V-Dem cont.", "Polity bin.", "JME lag", "Democracy",
                               "GDPPC (log)", "GDP",
                               "Military Expenditures (log)",  "CINC",
                               "Alliances",
                               "Rivalries (count)",
                               "Time", "Time$^{2}$", "Time$^{3}$"),
          column.separate = c(2, 1, 1, 1,1),
          column.labels = c("Alt. Democracy", "JME Lag", "GDP", "CINC"),
          dep.var.caption="",
          dep.var.labels.include = FALSE,
          add.lines=list(c('Countries',161,161,161,161,161)),
          model.numbers = FALSE,
          coef=list(contDemCoef, polityDemCoef, jmeLagCoef, gdpCoef, cincCoef),
          se=list(contDemSE, polityDemSE, jmeLagSE, gdpSE, cincSE),
          p.auto=TRUE,
          digits=3,
          star.cutoffs = .05,
          model.names=FALSE,
          omit.stat = c("f", "ser", 'bic'),
          notes=c("$^{*}$p$<$0.05; SEs clustered on 161 countries; all predictors lagged one year."),
          notes.append=FALSE,
          notes.align="l",
          table.placement = "H")

