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
#Use the updated data



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




##Preview the data




##Plot the data




##fit a linear model using the post as the exposure of interest




#Potential problems 


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




##Preview the data




##Plot the data




##fit a linear model using the post as the exposure of interest




#Potential problems 


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




##Preview the data




##Plot the data




##fit a linear model using the post as the exposure of interest




#Potential problems 


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







#---Step 2: Implementing the DID analysis method-----
# Fit a linear model with the treated, post and their interaction
# Don't forget to add any other time-varying covariates
# What is the regression coefficient of the interaction term

