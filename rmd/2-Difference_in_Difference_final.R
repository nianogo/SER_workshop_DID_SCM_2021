#------------------------------------------------------------------------------#
#--------An introduction to Difference-in-difference and Synthetic control------ 
#-------------------------methods for Epidemiologists"--------------------------
#-----------"Society for Epidemiologic Research (SER) Workshop: Part 2/3"------#
#-----------------------------Date:07/09/21------------------------------------#
#Roch Nianogo (niaroch@ucla.edu) & Tarik Benmarhnia (tbenmarhnia@health.ucsd.edu)
#------------------------------------------------------------------------------#

#------------------------------------Checking the directory---------------------
getwd()


#-------------------------------------Installing packages-----------------------


if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
} # a nice package to load several packages simultaneously


p_load("tidyverse","magrittr","broom",        #manipulate data
       "mise",                                #directory managment
       "here",                                #mise en place
       "Synth", "gsynth",                     #synthetic control
       "panelView", "lme4", "estimatr",       #multi-level model
       "ggdag",                               #draw a Directed Acylic Diagram (DAG)
       "gtsummary")                           #for tables

#---------------------------------Mise en place---------------------------------
mise()
#remove(list=ls())



#---------------------------------Loading the data------------------------------

#First load the data
mydata <- read_csv(here("data", "sim_data.csv"))
year_policy <- 2000

mydata <- mydata %>% 
  mutate(year_rec = year - year_policy,
         post     = ifelse(year>=year_policy,1,0),
         treated  = ifelse(state %in% c("Alabama",  "Alaska", 
                                        "Arizona", "Arkansas", "California"), 1,0),
         treatedpost = treated*post)


#-----------------------------------------Analysis------------------------------

#-----------------------#
#----Pre-post designs----
#-----------------------#
#Pre-post analysis (before-after) no control group: One state, two time points

#Subset the data to California and the years 1995 and 2005
#Preview the data
#Plot the data
#Fit a linear model to estimate the effect of the policy on the outcome y
#What are potential problems

##Create the data
dt <- mydata %>%
  filter(state=="California",
         year %in% c(1995, 2005)) 


##Preview the data
head(dt)


##Plot the data
dt %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

##fit a linear model using the post as the exposure of interest
fit <- lm(y ~ post + xi + xt + xit, data=dt)
tidy(fit)

# Potential problems
#-------------------#
# i.	The years before and after the policy have been implemented are picked arbitrarily
# ii.	It is impossible to get standard errors
# iii.	It is impossible to adjust for covariates
# iv.	This model would be wrong if the parallel trend assumption is violated
# v.	This model would be wrong if the common assumption is violated
# Conclusion 
# As expected the estimate is biased and standard errors inestimable


#----------------------------------#
#----Controlled Pre-post designs----
#----------------------------------#
#2 Pre-post analysis (before-after) with a control group: Two states, two time points


#Subset the data to California and Georgia and the years 1995 and 2005
#Preview the data
#Plot the data
#Fit a linear model to estimate the effect of the policy on the outcome y
#What are potential problems

##Create the data
dt1 <- mydata %>% 
  filter(state=="California" | state=="Georgia",
         year %in% c(1995, 2005))


##Preview the data
head(dt1)


##Plot the data
dt1 %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

##fit a linear model using the post as the exposure of interest
fit <- lm(y ~ treated*post, data=dt1) 


# Potential problems 
#--------------------#
# i.	The years before and after the policy have been implemented are picked 
# arbitrarily 
# ii.	The control units are picked arbitrarily
# iii.	It is impossible to get standard errors
# iv.	It is impossible to adjust for covariates
# v.	This model would be wrong if the parallel trend assumption is violated
# vi.	This model would be wrong if the common assumption is violated

# Conclusion 
# As expected the estimate is biased and standard errors inestimable



#---------------------------------------#
#----Interrupted Time Series Designs----
#---------------------------------------#
###Interupted time series: One state, Multiple time points


#Subset the data to California
#Preview the data
#Plot the data
#Fit a linear model to estimate the effect of the policy on the outcome y
#What are potential problems

##Create the data
dt2 <- mydata %>% 
  filter(state=="California")


##Preview the data
head(dt2)


##Plot the data
dt2 %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

##Fit the model
fit <- glm(y ~ year_rec + post + post*year_rec + xi + xt + xit, data=dt2)
tidy(fit)

# Potential problem (s)
#----------------------#
# This model would be wrong if the common assumption is violated 
# (that is if an event occurs in the treated unit at or after the time of 
# the policy that is not attributable to the policy)
# Conclusion 
# As expected the estimate is slightly biased 

#-------------------------------------------------#
#-----Controlled Interrupted Time Series (CITS)----
#----or Difference-in-Difference Designs-----------
#-------------------------------------------------#

#---Step 1: Checking the parallel trends assumption-----

#Check the parallel trends assumptions
#To do so
#restrict the data to the period before the policy
#Run a linear regression of the outcome on time-varying
#covariates and on an interaction term between treated indicator and year


p_load("estimatr") #for the lm_robust() function
pretrend_data <- mydata %>% 
  filter(post == 0)

res_pretrend <-  lm(y ~ treated*year + xit + xt + xi, 
                    data = pretrend_data)
summary(res_pretrend)


#the lm_robust procedure is best because SE are correctly estimated
res_pretrend <- lm_robust(y ~ treated*year + xit, data = pretrend_data,
                          fixed_effects=state,
                          clusters = state, se_type = "stata")
summary(res_pretrend)


#---Step 2: Implementing the DID analysis method-----
# Fit a linear model with the treated, post and their interaction
# Don't forget to add any other time-varying covariates
# What is the regression coefficient of the interaction term

#Method 1: Lm_robust()
p_load("estimatr")
all <- lm_robust(y ~ treated:post + post + xt + xit, 
                 data = mydata, fixed_effects=state,
                 clusters = state, se_type = "stata")
all


#Method 2: lmer(): Multilevel
fit <- lmerTest::lmer(y ~ treated + post + treated:post + xt + xi + xit + (1| state) + (1 | year), 
                      data = mydata, 
                      REML = T)
summary(fit)

#Other specifications
fit2 <-  lm(y ~  xit + xi + xt + factor(state) + factor(year) + treated:post, data = mydata)
summary(fit2)

p_load("lme4", "lmerTest")
fit3 <-  lmerTest::lmer(y ~ treated:post + xit + xi + xt + (1| state) + (1 | year), 
                        data = mydata, REML = T)
summary(fit3)

p_load("geepack")
res_model <- geeglm(y ~ treated*post + xi + xt + xit, 
                    family=gaussian("identity"), 
                    data=mydata %>% mutate(id=state_num), 
                    id=id)
tidy(res_model) #SE are a bit larger