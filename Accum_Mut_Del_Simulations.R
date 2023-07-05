library(IBMPopSim)
library(ggplot2)
library(cowplot)
library(tidyr)
library("viridis")
library(dplyr)
library(latex2exp)
library(plyr)


#We begin by creating the model with the package IBMPopSim

### Parameters
s0 = 0.9
s = 0.1
h = 0.25
shom = s
shet = h*s
mu_per_bp = 10^(-8)
mu = 0.1
epsilon = 0.1
birth_rate = 1
r=0.05

#With C++
params <- list("birth_rate" = birth_rate, "s0" = s0, "shom" = shom, "shet" = shet, "mu" = mu, 
               "epsilon" = epsilon, "rec" = r)



### Model creation

# Events definition, based on the possible events described in the main text

death_event <- mk_event_individual(type = "death", 
                                   intensity_code = "result = s0 + shom * I.w + shet*(I.xa+I.xb-2*I.w) ;")


birth_event <- mk_event_individual(type = "birth", intensity_code = "result = birth_rate;")

birth_rec_kernel_code <- "int a = CUnifInt(0,3);
double theta = CUnif(0,1);
if (CBern(rec)==1) {
if (a==1) {newI.xb = theta*I.xb + (1-theta)*I.xa;
newI.w=theta*I.w + (1-theta)*I.xa;}
else if (a==2) {newI.xa = theta*I.xa + (1-theta)*I.xb;
newI.w = theta*I.w + (1-theta)*I.xb;}
else if (a==3) {newI.xa = theta*I.xa + (1-theta)*I.xb;
newI.xb = theta*I.xb + (1-theta)*I.xa;}
}
"

birth_event_rec <- mk_event_individual(type = "birth", intensity_code = "result = birth_rate;", 
                                       kernel_code = birth_rec_kernel_code)


mutation_event <- mk_event_individual(type = "swap", intensity_code = "result = mu;", 
                                      kernel_code = "if (CBern(0.5) == 0) {I.xa = I.xa + epsilon*(1-I.xa);I.xb = I.xb;I.w = I.w + epsilon*(I.xb - I.w);} else {I.xa = I.xa;I.xb = I.xb+ epsilon*(1-I.xb);I.w = I.w + epsilon*(I.xa - I.w);}")

## Creation of the initial population (needed for the model creation)

#Initial size
init_size <- 1000
#Threshold to set initial values for x_a and x_b. Uniform distribution in [0, equ]
equ = 1

#Definition of the population, with three characteristics
xa = runif(init_size, min = 0, max = equ)
xb = runif(init_size, min = 0, max = equ)
#The third characteritic must satisfy max(0, x_a+x_b-1) < w < min(x_a, x_b)
w = rep(0, init_size)
for (i in 1:init_size){
  x = xa[i]
  y = xb[i]
  borne_inf = max(0, x+y-1)
  borne_sup = min(x, y)
  w[i] = runif(1, min = borne_inf, max = borne_sup)
}


#Creation of the initial pop
pop <- data.frame(birth = rep(0, init_size), death = NA, xa , xb, w)


# Model computation

model_nonrecombining <- mk_model(characteristics = get_characteristics(pop),
                                 events = list(death_event, birth_event, mutation_event), 
                                 parameters = params)


model_recombining <- mk_model(characteristics = get_characteristics(pop),
                              events = list(death_event, birth_event_rec, mutation_event), 
                              parameters = params)


### Settings for population simulation

#Specify bounds on the rates. Sarah Kaakai said that there was no need to have a precise bound
#In this program, we don't change the parameters. But it is possible to lauch simulations for different parameters
# without re-computing the model. Only need to lauch popsim with a different parameter vector specified
my_events_bounds <- c("death" = params$s0 + params$shet + params$shom, "birth" = params$birth_rate, "swap" = params$mu)


#===================
#Accessory functions
#===================


#Function that extracts the population size at given times
#times should be a vector containing the times at which we wish to obtain the pop size
#Data_total is supposed to be the data returned by IBMPopSim with a maximum time (not a vector of times)
population_sizes <- function(Data_total, times){
  sizes = c(1000)
  for (t in times){
    P = population_alive(Data_total$population, t)
    s = length(P$birth)
    sizes = c(sizes, s)
  }
  return(sizes)
}



#Function that plots the distribution of individuals traits (heterozygosity and homozygosity)
#This function is used to plot the graphs used in figure 3.4 of the manuscript.
#Data_total : the total result of an IBMPopSim simulation, launched with a maximum time that must be 
#the last entry of the vector "times".
#times : a vector that contains the times at which the population size is known and has to be plotted
#Sizes_princ and sizes_sec contains the populations sizes. Sizes_princ contains what we want to plot,
#Sizes_sec is just there for comparison and to ensure that the scales are the same for all graphs of 
#figure 3.4
heatmap_plot_scaled <- function(Data_total, times, sizes_princ, sizes_sec, k){
  
  #Extracting the population alive at the kth time of 'times'
  P = population_alive(Data_total$population, times[k])
  df = data.frame(omega = P$w, xhet = P$xa + P$xb - 2*P$w)
  
  #And constructing the heatmap
  heatmap <- ggplot(df, aes(xhet, omega))  + 
    scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    theme_classic() +
    stat_bin2d(bins = 100) +
    scale_fill_viridis(option='turbo', direction = 1) +
    ggtitle(paste('Population at time ', times[k])) +
    xlab(TeX(r'($x_{het}$)')) + 
    ylab(TeX(r'($\omega$)'))+
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), 
          text=element_text(size=17), plot.title = element_text(face = "plain", size=13), 
          axis.title.x = element_text(size=30), axis.title.y = element_text(size=30))
  
  
  #Creating the plot for the sizes
  df_size <- data.frame(time = c(0,times[1:k]), size = sizes_princ[1:(k+1)])
  max_size = max(sizes_princ[k+1], sizes_sec[k+1])
  print(max_size)
  
  pop_size <- ggplot(df_size, aes(x=time, y=size)) + geom_point() + geom_line()+
    theme_classic()+ ylim(0,max_size)+
    xlab('Time')+
    ylab('Population size')+
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), 
          text=element_text(size=17), plot.title = element_text(face = "plain", size=13))
  
  p = ggdraw()+
    draw_plot(pop_size, scale = 0.4, y= 0.2, x=0.15)+
    draw_plot(heatmap)+
    theme(plot.background = element_rect(fill='white'))
  
  print('heatmap OK')
  #ggsave(paste("./Evol_Rec/HeatSize_Timemax_",time_max, "_Time_", times[k], ".png",sep=""))
  
  return(p)
}


#==========================
#To launch the simulations and obtain figures presented in chapter 3 of the manuscript
#==========================

#Choose a maximum time until which the simulation will run. /!\ If the population size growths 
#exponentially, the computing limit may be reached.
time_max=80
#Vector of times at which the distribution will be plotted
time_vector = c(1, seq(10, time_max, 10))

#To launch a trajectory, recombining and non-recombining
sim_out_rec <- popsim(model=model_recombining, population = pop, events_bounds = my_events_bounds, parameters = params, time = time_max, multithreading = TRUE)
A = sim_out_rec
sim_out_norec <- popsim(model=model_nonrecombining, population = pop, events_bounds = my_events_bounds, parameters = params, time = time_max, multithreading = TRUE)
B = sim_out_norec
#Compute the population size to plot it
A_sizes = population_sizes(A, time_vector)
B_sizes = population_sizes(B, time_vector)
print(A_sizes)
print(B_sizes)
#Heatmap of omega/xhet and evolution of size.
#Change the last digit to chose the index of time_vector you want to plot
heatmap_plot_scaled(A, time_vector, A_sizes, B_sizes, length(time_vector))
heatmap_plot_scaled(B, time_vector, B_sizes, A_sizes, length(time_vector))
#Saving last figure
ggsave(paste("./timemax_",time_max, "_PopRec_", time_vector[9], ".png",sep=""))


#===========================================
#To plot histogram of death rates
A_alive = population_alive(A$population, time_max)
nb_ind_rec = length(A_alive$birth)
D_rec = data.frame(Inv = rep(1, nb_ind_rec), d = A$arguments$parameters$s0 + A$arguments$parameters$shom*A_alive$w + A$arguments$parameters$shet*(A_alive$xa + A_alive$xb - 2*A_alive$w ))

B_alive = population_alive(B$population, time_max)
nb_ind_norec = length(B_alive$birth)
D_norec = data.frame(Inv = rep(0, nb_ind_norec), d = B$arguments$parameters$s0 + B$arguments$parameters$shom*B_alive$w + B$arguments$parameters$shet*(B_alive$xa + B_alive$xb - 2*B_alive$w ))


df_d = rbind(D_rec, D_norec)
df_d$Inv <- factor(df_d$Inv, levels = c(1,0))
mu <- ddply(df_d, "Inv", summarise, grp.mean=mean(d))

ggplot(df_d, aes(x=d, fill = Inv)) + 
  geom_histogram(position="identity", bins=60, alpha = 0.6) +
  scale_fill_manual(values = c( "#ffb53e", "#5b06ac"))+
  labs(title=paste("Histogram of death rates for time", time_max), x="Death rate", y="Count", fill="R")+
  theme_classic()+
  theme(text=element_text(size=13), plot.title = element_text(face = "plain", size=13))+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=as.factor(Inv)), linewidth = 1.2, show.legend=FALSE)+
  scale_color_manual(values = c( "#ffb700", "#3b0073"))

ggsave("test.png", width=20, units="cm")


#============================================================ 

