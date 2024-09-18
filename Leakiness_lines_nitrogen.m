
% Code to generate supplement Figure 8

% Want growth rate (gNT) for varying leakiness in even environmnts. 
% Distinct panels for low and high reward rate. 


clear all
figure
clf

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of nonstructural carbon to mychorrhizal carbon pool
s = 0.01; %senesence of nonstructural carbon
lambda = 0.005; %loss of mycorrhizal carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.005; %fungus 1 mortality
m2 = 0.005; %fungus 2 mortality
sN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

u_bar = .5; %mean uptake of nitrogen by fungus
difference_val = .5; %severity of environment, high
u1_A = u_bar + difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = u_bar - difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = u_bar - difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = u_bar + difference_val;%uptake of Nitrogen by fungus 2 in environment type B

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

%density dependent mortality paramters 
d1_1 = .1;
d2_1 = .05;
d2_2 = .1;
d1_2 = .05;

leakiness_vals = [0:.025:1];

colorval = cmocean('deep', 6); %need this to get correct colors https://matplotlib.org/cmocean/
rtot_vals = [.2 .8]; 
env_periods = [365 365*2  365*4 365*6];

for rt = 1:2 %each panel is a diffferent reward rate 
    rtot = rtot_vals(rt) ;

    for k = 1:length(env_periods); %each panel has a few lines for different periods
        env_period = env_periods(k);
        propA = 0.5;
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;
        
        results = nan(1, length(leakiness_vals));
        for i = 1:length(leakiness_vals)
            leakiness = leakiness_vals(i);

            tspan = [1 80000];
            if tspan(2)-env_period*4 < 8000
               tspan = [1 env_period*8];
            end

            %simulate
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, lambda, rtot*(1-leakiness), rtot*leakiness, rtot*leakiness, rtot*(1-leakiness), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            biomass_val = mean(final_res(1,:)); %store average photosynthetic biomass 
            nitrogen_val = mean(final_res(5,:)); %store average tree nitrogen  

            results(i) = nitrogen_val ;

        end
        %plot 
        subplot(1,2,rt)
        hold on
        plot(leakiness_vals, g.*results, 'linewidth', 2, 'color', colorval(k+1, :));

    end
end

%some titles and labels 
subplot(1,2,1)
title('Low reward rate')
ylim([10 18])
ylabel('Tree carbon growth rate (gN_T)')
xlabel('Leakiness')
box on

legend({'1 year'; '2 years'; '4 years'; '6 years'}, 'Location', 'southwest')
subplot(1,2,2)
ylim([10 18])
xlabel('Leakiness')
title('High reward rate')
box on
