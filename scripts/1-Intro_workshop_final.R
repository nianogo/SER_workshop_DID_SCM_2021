#------------------------------------------------------------------------------#
#--------An introduction to Difference-in-difference and Synthetic control------ 
#-------------------------methods for Epidemiologists"--------------------------
#-----------"Society for Epidemiologic Research (SER) Workshop: Part 1/3"------#
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

# This is a simulated data where a fictitious policy ban was implemented
# The policy was enacted in five states: Alabama,  Alaska, Arizona, Arkansas, California
# The policy was enacted in 2000
# The unit of analysis is the state
# y  = is the outcome
# xi = time-invariant variable (but varies across states)
# xt = time-varying variable (but is constant across states)
# xit = time-varying and unit-variable (i.e. varies by year and state)




#First load the data

mydata <- read_csv(here("data", "raw_data", "sim_data.csv"))


#Create some useful variables within the datasets:
#an indicator for after the policy has been implemented: called this "post"
#an indicator for states that have implemented the policy: called this "treated"
#an interaction term between posst and treated: called this treatedpost
#create a recentered year variable such that year_rec = 0 at the time of the policy

year_policy <- 2000

mydata <- mydata %>% 
  mutate(year_rec = year - year_policy,
         post     = ifelse(year>=year_policy,1,0),
         treated  = ifelse(state %in% c("Alabama",  "Alaska", 
                                        "Arizona", "Arkansas", "California"), 1,0),
         treatedpost = treated*post)

#-------------------------Explore the data structure----------------------------

head(mydata, 10)

glimpse(mydata)

summary(mydata)

p_load("skimr")
skim(mydata)

#-----------------------------------Plot a DAG------------------------------
p_load("ggdag")

coords <- list(
  y = c( 
    Post = 1, 
    Year = 0.5,
    Treat = 1.5,
    Treated_post =1,
    States = 2,
    C = 1.5,   
    Y = 1
  ),
  x = c( 
    Post = 1, 
    Year  = 2,
    Treat = 2,
    Treated_post =3,
    States = 3,
    C = 4,   
    Y = 5
  )
)



dag0 <- dagify(Y ~ Treated_post + States + Year,
               Post  ~ Year,
               Treat ~ States,
               Treated_post ~ Post + Treat,
               coords=coords) 

ggdag(dag0, node_size = 18, node = T) + 
  geom_dag_point(color="white") +
  geom_dag_text(color="black", size=5) +
  theme_dag()


#----------------------------Visualize the data and policy----------------------

#--------------------------#
#----Treatment overview-----
#--------------------------#

mydata <- mydata %>% 
  mutate(year_rec = as.integer(year_rec)) %>%
  as.data.frame() # need to convert to a data.frame for some function to work


#----Using PanelView---
p_load("panelView")
panelView(y ~ treatedpost, data = mydata,
          index = c("state","year_rec"), 
          pre.post = TRUE) 


#----Using Panel Match---
p_load("PanelMatch")
DisplayTreatment(unit.id = "state_num",
                 time.id = "year_rec", legend.position = "none",
                 xlab = "year", ylab = "state",
                 treatment = "treatedpost", data = mydata)


#----Map of US with place where is implemented---
p_load("maps", "ggthemes")
us_states <- map_data("state")

mydata2014 <- mydata %>% 
  filter(year==2000) %>% 
  mutate(region = str_to_lower(state))

mydata_maps <- left_join(mydata2014, us_states, by = "region")

p <- ggplot(data = mydata_maps,
            aes(x = long, y = lat, 
                group = group, fill = factor(treated))) +
  labs(title = "Smoking bans in the US in 2000", fill = NULL) +
  geom_polygon(color = "gray90", size = 0.1) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)  +
  theme_map()  +
  guides(fill = guide_legend(nrow = 3)) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")
p
#Note that Alaska (which is treated) is not shown from the map

#-------------------------------------------------#
#----Sample tables and Covariate balance plots----
#--------------------------------------------------#

##Table 1
p_load("gtsummary")
n_treated0 <- mydata %>% 
  filter(treated == 0) %>% 
  select(state) %>% 
  n_distinct()

n_treated1 <- mydata %>% 
  filter(treated == 1) %>% 
  select(state) %>% 
  n_distinct()

tab1 <- mydata %>%
  filter(post==0) %>% 
  select(c("xit","xt", "xi", "treated","y")) %>% 
  mutate(treated= case_when(treated==1~"Treated",
                            TRUE~"Control")) %>% 
  tbl_summary(
    missing = "no",
    by ="treated",
    type = list(everything() ~ "continuous"),
    digits = list(everything() ~ 2),
    statistic = list(everything()~"{mean}")
  ) %>% 
  modify_spanning_header(c("stat_1", "stat_2") ~ "Table 1") %>% 
  modify_header(update = list(
    label ~ '**Characteristic**',
    stat_1 ~ '**Control**, N = {n_treated0}',
    stat_2 ~ '**Treated**, N = {n_treated1}'
  )) 
tab1 %>% as_flex_table()


##Covariate balance plot
##---By treated group
p1 <- mydata %>%
  filter(post==0) %>% 
  select(xt, xi, xit, y, treated) %>% 
  mutate(treated= case_when(treated==1~"Treated",
                            TRUE~"Control")) %>% 
  group_by(treated) %>% 
  group_modify(~ {.x %>% map_dfr(mean)}) %>% 
  pivot_longer(cols = -treated,
               names_to = c("variable"),
               values_to = "mean") %>% 
  ggplot(aes(x=treated, y=mean, fill=factor(treated))) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label = round(mean,2)), 
            position = position_stack(0.5), 
            size=3, 
            color = "black")+
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Checking for imbalance in variables pre-policy",
       y = "Mean",
       x = "Variables",
       fill = "Treatment status")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 
p1
#notice xt is constant across units

##---By post period
p2 <- mydata %>%
  filter(treated==0) %>% 
  select(xt, xi, xit, y, post) %>% 
  mutate(post= case_when(post==1~"After",
                         TRUE~"Before")) %>% 
  group_by(post) %>% 
  group_modify(~ {.x %>% map_dfr(mean)}) %>% 
  pivot_longer(cols = -post,
               names_to = c("variable"),
               values_to = "mean") %>% 
  ggplot(aes(x=post, y=mean, fill=factor(post))) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label = round(mean,2)), 
            position = position_stack(0.5), 
            size=3, 
            color = "black")+
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Checking for imbalance in variables in control units",
       y = "Mean",
       x = "Variables",
       fill = "Time")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 
p2
#notice xi is constant across time

#------------------------#
#----Time trends plots----
#------------------------#


#For all the data
mydata %>% 
  ggplot(aes(x=year, y=y, group=state)) + 
  annotate("rect", fill = "gray", alpha = 0.5,
           xmin = 2000, xmax = 2010,
           ymin = -Inf, ymax = Inf) +
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       color = "Treatment") +
  geom_line(aes(color=factor(treated)), size=0.5) +
  scale_color_discrete(labels=c("Control", "Treated")) +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 


#For two states
mydata %>% 
  filter(state %in% c("California", "Georgia")) %>% 
  ggplot(aes(x=year, y=y, group=state)) + 
  annotate("rect", fill = "gray", alpha = 0.5,
           xmin = 2000, xmax = 2010,
           ymin = -Inf, ymax = Inf) +
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       color = "Treatment") +
  geom_line(aes(color=factor(treated)), size=0.5) +
  scale_color_discrete(labels=c("Control", "Treated")) +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 


#On average
mydata %>% 
  group_by(year, treated) %>% 
  summarise(y=mean(y),.groups="keep") %>% 
  ggplot(aes(x=year, y=y, group=treated, color = factor(treated))) + 
  annotate("rect", fill = "gray", alpha = 0.5,
           xmin = 2000, xmax = 2010,
           ymin = -Inf, ymax = Inf) +
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  scale_color_discrete(labels=c("Controls", "Treated")) +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 
