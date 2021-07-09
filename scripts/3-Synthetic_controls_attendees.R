#------------------------------------------------------------------------------#
#--------An introduction to Difference-in-difference and Synthetic control------ 
#-------------------------methods for Epidemiologists"--------------------------
#-----------"Society for Epidemiologic Research (SER) Workshop: Part 3/3"------#
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



#-----------------------------------------Analysis------------------------------

#----------------------------------------------------------------#
#-------Traditional Synthetic Control Method (Abadie)-------------
#----------------------------------------------------------------#

#the data needs to be a data.frame (as opposed to a tibble) 
#to be able to be used in the Synth package.
#One can use the as.data.frame() command

mydata <- mydata %>% 
  mutate(year_rec = as.integer(year_rec)) %>%
  as.data.frame() # need to convert to a data.frame for some function to work


# create a list of treated states
list_int <- mydata %>% 
  filter(treated==1) %>% 
  select(state_num) %>% 
  distinct() %>% 
  as.matrix() %>% 
  as.vector()
list_int

# create a list of control states
list_control <- mydata %>% 
  filter(treated==0) %>% 
  select(state_num) %>% 
  distinct() %>% 
  as.matrix() %>% 
  as.vector()
list_control

# create some handy object to be used in the Synth package
year_start  = 1961
year_end    = 2010
year_policy = 2000
year_policy_prior = year_policy - 1


p_load("Synth")
#----Step 1: Prep the data
dataprep.out <-
  dataprep(foo = mydata,
           predictors = c("xit"), #include time-varying and unit-varying variabless
           predictors.op = "mean",
           dependent     = "y",
           unit.variable = "state_num",
           time.variable = "year",
           unit.names.variable = "state",
           treatment.identifier  = 1, #cannot be several states
           controls.identifier   = list_control,
           time.predictors.prior = c(year_start:year_policy_prior),
           time.optimize.ssr     = c(year_start:year_policy_prior),
           time.plot             = c(year_start:year_end))


#----Step 2: Run the Synth command
synth.out <- synth(dataprep.out)

#----Step 3: Get result tables
print(synth.tables   <- synth.tab(
  dataprep.res = dataprep.out,
  synth.res    = synth.out)
)

# see thee weights
round(synth.out$solution.w,2)


#----Step 4: Obtain the estimate using a DID
##--Estimate the average (weighted) outcome in the synthetic control prior to the policy

y_ctrl_pre <- dataprep.out$Y0plot%*%synth.out$solution.w %>% 
  as.data.frame() %>% 
  rownames_to_column(var="year") %>% 
  rename(control=w.weight) %>% 
  filter(year %in% c(year_start:year_policy_prior)) %>% 
  group_by() %>% 
  summarise(mean=mean(control)) %>% 
  as.matrix() %>% 
  as.vector()
y_ctrl_pre

##--Estimate the average (weighted) outcome in the synthetic control after the policy
y_ctrl_post <- dataprep.out$Y0plot%*%synth.out$solution.w %>% 
  as.data.frame() %>% 
  rownames_to_column(var="year") %>% 
  rename(control=w.weight) %>% 
  filter(year %in% c(year_policy:year_end)) %>% 
  group_by() %>% 
  summarise(mean=mean(control)) %>% 
  as.matrix() %>% 
  as.vector()
y_ctrl_post

##--Estimate the average outcome in the treated state prior to the policy
y_treat_pre <- dataprep.out$Y1plot %>% 
  as.data.frame() %>% 
  rownames_to_column(var="year") %>% 
  rename(treated=`1`) %>% 
  filter(year %in% c(year_start:year_policy_prior)) %>% 
  group_by() %>% 
  summarise(mean=mean(treated)) %>% 
  as.matrix() %>% 
  as.vector()
y_treat_pre

##--Estimate the average outcome in the treated state after the policy
y_treat_post <- dataprep.out$Y1plot %>% 
  as.data.frame() %>% 
  rownames_to_column(var="year") %>% 
  rename(treated=`1`) %>% 
  filter(year %in% c(year_policy:year_end)) %>% 
  group_by() %>% 
  summarise(mean=mean(treated)) %>% 
  as.matrix() %>% 
  as.vector()
y_treat_post 


##--Estimate a simple difference-in-difference
diff_post = y_treat_post - y_ctrl_post
diff_pre = y_treat_pre - y_ctrl_pre
diff_pre

DID = diff_post - diff_pre
DID


#----Step 5: Plot the results
##--Plot treated and synthetic state
path.plot(synth.res    = synth.out,
          dataprep.res = dataprep.out,
          Ylab         = c("Y"),
          Xlab         = c("Year"),
          Legend.position = c("topleft"))

abline(v   = year_policy,
       lty = 2)


##--Plot gap between treated and synthetic state
gaps.plot(synth.res    = synth.out,
          dataprep.res = dataprep.out,
          Ylab         = c("Gap"),
          Xlab         = c("Year"),
          Main         = ""
)
abline(v   = year_policy,
       lty = 2)


#------------------------------------------------------#
#------Preliminary notions (regression modeling--------
#----through matrix multiplication and optimization----
#------------------------------------------------------#

#simulate some data
set.seed(1)
x <- matrix(rbinom(1000, 1, prob=0.3), nrow=100, ncol=1, byrow=T)
y <- matrix(rnorm(1000, 2 + 3*x, sd=1), nrow=100, ncol=1, byrow=F)
int <- rep(1,length(y))
xint <- cbind(int,x)
xint

#Method 1: Ordinary Least Squares (OLS) through simple linear regression
mod1 <- lm(y~x)
tidy(mod1)

#Method 2: Ordinary Least Squares (OLS) through matrix multiplication
solve(t(xint) %*% xint) %*% t(xint) %*% y

#Method 3: Ordinary Least Squares (OLS) through optimzation
p_load("CVXR")
n <- ncol(xint)
beta <- Variable(n)
objective <- Minimize(p_norm((y - xint%*%beta), p=2))
prob <- Problem(objective)
CVXR_result <- solve(prob)
cvxrBeta <- CVXR_result$getValue(beta) 
cvxrBeta

#------------------------------------------------------#
#--------------Manual Synthetic Control Method---------
#------------------------------------------------------#

#---Step 1: 
# Aggregate the data by state, treated and post
# Transpose the data to obtain a unique row per state
# Create difference variable y_diff = y_1 - y_0
# Obtain the weights through optimization (steps shown below)
# Apply the weights to check the covariate balance
# fit a weighted regression model





#---Step 2: Obtain the weights through optimization
p_load("CVXR") #optimization package: convex optimzation

##--Obtaining the value of the xit variable for treated and control
cov_treat <- df00 %>% 
  ungroup() %>% 
  filter(state=="Alabama") %>% #to do one treated state at a time
  select(xit_0) %>% 
  as.matrix() %>% 
  t()
cov_treat

dim(cov_treat)
cov_treat %>% t()

cov_ctrl <- df00 %>% 
  ungroup() %>% 
  filter(treated==0) %>% 
  select(xit_0) %>% 
  as.matrix() %>% 
  t()

cov_ctrl

dim(cov_ctrl)
cov_ctrl %>% t()


##--obtaining the weights through optimization
# The SCM is formed by finding the vector of weights W
# that minimizes the imbalance between the treated unit and 
# weighted average of the controls across a set of variables
# C (e.g. pre-intervention covariates)
# Essentially in a SCM we try to find W such that 
# (C1 - C0W)'V(C1- C0W) is minimized
#W being positive and summing to 1 to avoid extrapolating (Abadie et al 2010). 

v <- Variable(c(ncol(cov_ctrl), ncol(cov_treat)))
obj <- p_norm((cov_treat - cov_ctrl%*%v), p=2)
objective <- Minimize(obj)
prob <- Problem(objective, list(v>=0, sum(v)==1)) #we add in the constraints
CVXR_result <- solve(prob)
synth_weights <- CVXR_result$getValue(v)
round(synth_weights,3)

##--comparing the weights with those of the synthetic control package
cbind(round(synth_weights,3), synth.tables$tab.w)


##--Apply the weights to check covariate






##--Implement a DID





##--Fit a weighted regression model of y_diff on the treated indicator
myweights <- c(1, synth_weights) # weight of 1 for the treated group
myweights 

df00a <- df00 %>% 
  filter(state=="Alabama" | treated == 0) %>% 
  mutate(weights = myweights) 

head(df00a,20)



#---------------------------------------#
#----Generalized Synthetic Control-------
#---------------------------------------#
#https://yiqingxu.org/software/gsynth/gsynth_examples.html
p_load("gsynth")
y <- gsynth(y ~ treatedpost + xit, 
            data = mydata,  EM = F, index = c("state","year"), 
            inference = "parametric", se = TRUE,
            nboots = 2,  r = c(0, 5), CV = TRUE, force = "two-way", parallel = FALSE)
y1 <- round(data.frame(y$est.avg),2)
y1


plot(y)
plot(y, type = "counterfactual", raw = "none", main="")
plot(y, type = "counterfactual", raw = "all")
plot(y, type = "ct", raw = "none", main = "", shade.post = FALSE)

