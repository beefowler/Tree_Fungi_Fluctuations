This repository has been produced in association with the manuscript, 
"Environmental fluctuations promote host reward strategies
that maintain partner diversity in multispecies mutualisms"
which was accepted to the American Naturalist August 30, 2024.
  
Authors: Bethany L.F. Stevens, Kristen Howard, Laura M. Bogar, Holly V. Moeller. 

contact bstevens@ucsb.edu for questions. 

Summary: In this work, we present a model of a simple-tree fungal mutualism in which partner quality changes over time. The model has two distinct environmental conditions and host that can employ carbon reward strategies with varying degrees of preference between higher and lower-quality fungus. We explore the dynamics of this system across environmental regimes and with different reward strategies in order to identify patterns in tree fitness and fungal community composition. 

All scripts in this repository were written and edited by BLFS based on earlier versions in R scripts and pseudocode laid out by all authors. 


Prerequisites: 
Some of the figures make use of the boundedline.m function by Kelly Kearney available here
https://github.com/kakearney/boundedline-pkg
Generated for the manuscript with version 1.4.01.

Some of the figures also rely on the cmocean colormaps by Kristen Thyng. https://matplotlib.org/cmocean/

Both of these packages are contained in subfolders in this repository for ease of use, in agreement with their associated licenses, but we had no role in creating those packages nor have we modified them specifically for our work. 

All files are run in MATLAB unless otherwise indicated. 
----------------------------------------

Experiments: 

leaky_or_loyal_coexistence.m 
	This is the function that contains the main system of differential equations laid forth in equations 1-5 in the manuscript (in "The Model" section just before "Environmental Forcing"). All inputs are explained in detail at the top of the script. Briefly, given parameter values and the state of the system at input time t, the function outputs the trajectory of the system as the derivative of each of the state variables at that time. 
 	This function is called by all subsequent scripts. 

Constant_reward_rates.m 
	This script allows us to investigate the behavior of the system under different environmental regimes when trees reward the two fungal partners at constant rates (in proportions of 1:0, 1:1, and 0:1). This script uses ode45 to simulate the system of differential equations under various values of alpha and delta (environmental parameters) for fixed biological parameter values. It stores the variance in the tree carbon pool in the final three periods of the simulation and plots the variance across environmental regimes as well as 3 example trajectories over time. 
	This script generates Figure 2 in the Manuscript. 

Leakiness_lines.m 
	This script allows us to compare the success of trees with "leaky" reward strategies of different degrees. Again, we use ode45 to simulate the system of differential equations under various environmental regimes (this time varying T and alpha only) for different values of the total reward rate and leakiness value l. The script stores only the average nonstructural tree carbon pool in the final 3 periods of each simulation for plotting. 
	This script generates Figure 3 in the Manuscript. 
 
Responsive_reward_strategies.m 
	This script provides more detail on the dynamics of the system when trees have a subset of leaky reward strategies. In particular, it calculates and plots the variation of the nonstructural tree carbon pool and fungal carbon pools when trees have values of l = 0, 0.25 and 0.75 across different environmental regimes (varying alpha).
	This script generates Figure 4 in the Manuscript.

Bifurcation_Supplement_Grouped.nb
	This script is run in Wolfram Mathematica.
	This script performs an equilibrium bifurcation analysis for the simplified model presented in part 1 of the supplement (Model dynamics in a constant environment). The system of equations is symbolically defined and then equilibria are calculated for a sweep of values for each parameter in the model. In each panel, all other parameter values are held constant. 
	This script generates Figures S1 in the Supplement.

ExtinctionStability_Supplement.m 
	This script explores the basin of attraction of the extinction equilibrium for the model system in a constant environment. In particular, we ask, when the system has very low initial values of fungus carbon, F(0), do the fungus and the tree go extinct or grow to their non-zero equilibrium? We hold all parameters constant except for reward rate and use ode45 to simulate the system with different initial conditions. 
	This script generates Figure S2 in the Supplement. 

Assymmetric_Fungus.m 
	This script addresses the question of whether our first result holds when fungal partners are not "mirror images" of each other, as explained in part 2 of the supplement, "Asymmetric demographic parameters". That is, we vary the mortality and efficiency parameters for one fungal partner while leaving the other constant. The script operates very similarly to constant_reward_rates.m in that for each set of biological parameters we want to explore, we simulate the system with ode45 for a range of possible values of alpha (0 to 1) and store the variation in the nonstructural carbon pool in the final three periods of the simulation. 
	This script generates Figure S3 in the Supplement. 

Simulations_supplement.m 
	We use this script to plot the values of all state variables in the system over time, since most of our other figures present carbon only. This script simulates the model using ode45 and plots the last few periods of dynamics for carbon pools, fungal biomass and nitrogen pools in three separate panels for each simulation. 
	This script generates Figure S4 in the Supplement. 

Partner_Preference_Supplement.m 
	This script demonstrates how the analysis presented in Figure 2 is relatively insensitive to the degree of preference between the two fungal partners. It operates along the same lines as constant_reward_rates.m, except with fungal reward ratios of 60:40, 53:48, and 51:49 instead of  1:0. 
	This script generates Figure S5 in the Supplement.

Initial_conditions_coexistence_supplement.m 
	This script caries out the analysis described in part 3 of the supplement, "Competition and sensitivity to initial conditions". It simulates the model across a range of environmental regimes (alpha = 0-1), for different combinations of competition parameters. In particular, we consider cases where interspecific competition is less than, equal to, and greater than intraspecific competition. We also use this script to test the dynamics of these simulations for two initial conditions, one in which the two fungal partners start at equal densities, and on where F1 has 10x as much carbon biomass as F2 at the start of the simulation. 
	This script generates Figure S6 in the Supplement. 

Leakiness_Contour_with_convergence.m 
	This script expands on the analysis presented in Leakiness_lines.m in that we here allow the environmental regime to vary in both evenness (alpha) and periods, and simulate the system when trees have a range of leakiness values. For each combination of parameters, the model is simulated with ode45 for 10000 days. At this point the script checks to see if the solution has converged, and if it has not, the simulation is continued for another 10000 days until the solution converges. We also check that at the end of the simulation both fungal carbon pools are above a minimum threshold of 0.01, to see if either went extinct. 
	This script generates Figure S7 in the Supplement. 

Leakiness_lines_nitrogen.m 
	This script tests whether our results are sensitive to our choice of the fitness proxy. In particular we wanted to see if the optimal leaky strategy would change if fitness were measured as the rate of carbon fixation rather than the pool size. We do this by modifying Leakiness_lines.m only to report the average value of the tree nitrogen pool rather than the nonstructural carbon pool (since photosynthetic growth rate is proportional to the size of this nitrogen pool). 
	This script generates Figure S8 in the Supplement. 

-------------------------------------------------------------------------------------------------


Data: 
Values for the data points plotted in each panel of each figure in the manuscript are included in the following csv files. These are the simulated values produced in the scripts above and described in detail in the figure captions. They are provided here for those who do not wish to run the model scripts to produce the simulated values. Column names and descriptions are included below. 


Figure2_TopPanelsData.csv 
Data included in upper 3 panels of Figure 2 only. Column 1 is x-coordinate. All curves can be drawn by plotting subsequent columns as y-coordinate against Column 1. 
	
alpha - value of parameter alpha, proportion of time spent in environment A, used as x-coordinate for all lines in these panels 

Panel_1_Bet_hedging - value of C_p (nonstructural carbon pool) in bet-hedging scenario in more severe environment (left panel). Black horizontal line in figure. 

Panel_1_F1_preference_25percentile - 25th percentile value of C_p (nonstructural carbon pool) when tree has F1 preference in more severe environment (left panel). Lower bound of pink shaded region. 

Panel_1_F1_preference_50percentile - Median value of C_p (nonstructural carbon pool) when tree has F1 preference in more severe environment (left panel). Thick pink curve. 

Panel_1_F1_preference_75percentile - 75th percentile value of C_p (nonstructural carbon pool) when tree has F1 preference in more severe environment (left panel). Upper bound of pink shaded region. 
	
Panel_1_F2_preference_25percentile - 25th percentile value of C_p (nonstructural carbon pool) when tree has F2 preference in more severe environment (left panel). Lower bound of pink shaded region. 	

Panel_1_F2_preference_50percentile - Median value of C_p (nonstructural carbon pool) when tree has F2 preference in more severe environment (left panel). Thick pink curve. 
	
Panel_1_F2_preference_75percentile - 75th percentile value of C_p (nonstructural carbon pool) when tree has F2 preference in more severe environment (left panel). Upper bound of pink shaded region. 
	
Panel_2_Bet_hedging - value of C_p (nonstructural carbon pool) in bet-hedging scenario in moderate environment (middle panel). Black horizontal line in figure. 

Panel_2_F1_preference_25percentile - 25th percentile value of C_p (nonstructural carbon pool) when tree has F1 preference in moderate environment (middle panel). Lower bound of pink shaded region. 

Panel_2_F1_preference_50percentile -  Median value of C_p (nonstructural carbon pool) when tree has F1 preference in moderate environment (middle panel). Thick pink curve. 
	
Panel_2_F1_preference_75percentile - 75th percentile value of C_p (nonstructural carbon pool) when tree has F1 preference in moderate environment (middle panel). Upper bound of pink shaded region. 
	
Panel_2_F2_preference_25percentile - 25th percentile value of C_p (nonstructural carbon pool) when tree has F2 preference in moderate environment (middle panel). Lower bound of pink shaded region. 

Panel_2_F2_preference_50percentile - Median value of C_p (nonstructural carbon pool) when tree has F2 preference in moderate environment (middle panel). Thick pink curve. 

Panel_2_F2_preference_75percentile - 75th percentile value of C_p (nonstructural carbon pool) when tree has F2 preference in moderate environment (middle panel). Upper bound of pink shaded region. 	

Panel_3_Bet_hedging - value of C_p (nonstructural carbon pool) in bet-hedging scenario in more mild environment (right panel). Black horizontal line in figure. 
	
Panel_3_F1_preference_25percentile - 25th percentile value of C_p (nonstructural carbon pool) when tree has F1 preference in more mild environment (right panel). Lower bound of pink shaded region. 

Panel_3_F1_preference_50percentile -  Median value of C_p (nonstructural carbon pool) when tree has F1 preference in more mild environment (right panel). Thick pink curve. 
	
Panel_3_F1_preference_75percentile - 75th percentile value of C_p (nonstructural carbon pool) when tree has F1 preference in more mild environment (right panel). Upper bound of pink shaded region. 
	
Panel_3_F2_preference_25percentile - 25th percentile value of C_p (nonstructural carbon pool) when tree has F2 preference in more mild environment (right panel). Lower bound of pink shaded region. 

Panel_3_F2_preference_50percentile - Median value of C_p (nonstructural carbon pool) when tree has F2 preference in more mild environment (right panel). Thick pink curve. 

Panel_3_F2_preference_75percentile - 75th percentile value of C_p (nonstructural carbon pool) when tree has F2 preference in more mild environment (right panel). Upper bound of pink shaded region. 	





Figure2_BottomPanelsData.csv 
Data included in lower 3 panels of Figure 2 only. Column 1 is x-coordinate. All curves can be drawn by plotting subsequent columns as y-coordinate against Column 1. These panels are simply time series showing the simulation values for nonstructural carbon, C_p, and Fungal carbon, F1 and F2. 

time - time in days since start of the simulation. 
Panel_1_C_p - Simulated value of nonstructural carbon variable in more severe environment (left panel)
Panel_1_F1 - Simulated value of Fungus 1 carbon pool in more severe environment (left panel)
Panel_1_F2 - Simulated value of Fungus 2 carbon pool in more severe environment (left panel)

Panel_2_C_p - Simulated value of nonstructural carbon variable  in moderate environment (middle panel)	
Panel_2_F1 - Simulated value of Fungus 1 carbon pool in moderate environment (middle panel)
Panel_2_F2 - Simulated value of Fungus 2 carbon pool in moderate environment (middle panel)	


Panel_3_C_p - Simulated value of nonstructural carbon variable in more mild environment (right panel)
Panel_3_F1 - Simulated value of Fungus 1 carbon pool in more mild environment (right panel)	
Panel_3_F2 - Simulated value of Fungus 2 carbon pool in more mild environment (right panel)





	
Figure3_Data.csv 
Data included in both panels of Figure 3. Column 1 is x-coordinate. All curves can be drawn by plotting subsequent columns as y-coordinate against Column 1.

leakiness - leakiness of host reward strategy. Proportion of reward sent to low quality fungal partner. used as x-coordinate for all lines in these panels 
Panel1_Cp_1yr - Average value of nonstructural carbon pool (Cp) from final 3 periods of simulations subject to environmental fluctuations with period 1 year and reward rate 0.2. 
Panel1_Cp_2yr - Average value of nonstructural carbon pool (Cp) from final 3 periods of simulations subject to environmental fluctuations with period 2 years and reward rate 0.2. 
Panel1_Cp_4yr - Average value of nonstructural carbon pool (Cp) from final 3 periods of simulations subject to environmental fluctuations with period 4 years and reward rate 0.2. 
Panel1_Cp_6yr - Average value of nonstructural carbon pool (Cp) from final 3 periods of simulations subject to environmental fluctuations with period 6 years and reward rate 0.2. 
Panel2_Cp_1yr - Average value of nonstructural carbon pool (Cp) from final 3 periods of simulations subject to environmental fluctuations with period 1 year and reward rate 0.8. 
Panel2_Cp_2yr - Average value of nonstructural carbon pool (Cp) from final 3 periods of simulations subject to environmental fluctuations with period 2 years and reward rate 0.8. 
Panel2_Cp_4yr - Average value of nonstructural carbon pool (Cp) from final 3 periods of simulations subject to environmental fluctuations with period 4 years and reward rate 0.8. 
Panel2_Cp_6yr - Average value of nonstructural carbon pool (Cp) from final 3 periods of simulations subject to environmental fluctuations with period 6 years and reward rate 0.8. 



Figure4_NoleakinessData.csv
Data included in leftmost panels of Figure 4 only. Column 1 is x-coordinate. All curves can be drawn by plotting subsequent columns as y-coordinate against Column 1.


alpha - value of parameter alpha, proportion of time spent in environment A, used as x-coordinate for all lines in these panels 

C_p_bet_hedging	
C_p_switching_25percentile	
C_p_switching_50percentile	
C_p_switching_75percentile	
Fungus_bet_hedging	
F1_switching_25percentile	
F1_switching_50percentile	
F1_switching_75percentile	
F2_switching_25percentile	
F2_switching_50percentile	
F2_switching_75percentile


 