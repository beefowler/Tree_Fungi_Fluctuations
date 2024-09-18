%Code to generate supplementary figure 6
% importance of density dependent mortality for avoiding sensitivity to initial
% conditions

clear all 
figure

% set parameter values randomly 
% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.004; %fungus 1 mortality
m2 = 0.004; %fungus 2 mortality
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

u_bar = .5; %mean uptake of nitrogen by fungus
difference_val = .5; %severity of environment, high
u1_A = u_bar + difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = u_bar - difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = u_bar - difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = u_bar + difference_val;%uptake of Nitrogen by fungus 2 in environment type B

rtot = .8; 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(4) = 1; %F2
x0(5) = 1; %N


% Set timespan and environment conditions during timespan 
tspan = [1 8000]; 
env_period = 365; 

%initialize results
propA_vals = 0:.02:1; 
u2_Bvals = 0:.1:1;
results = nan(3,length(propA_vals)); %mean value 
results2A = nan(3, length(propA_vals)); %median and iqr
results2B = results2A; 
results3 = results2B; 


preference=1; 

d1_1 = .1;
d2_2 = .1;

count = 1; %useful for plotting 
for init = 1:2 %two rows of different initial conditions 

    if init == 1
        x0(3) = 1; %F1
    else
        x0(3) = 10; %F1
        tspan = [1 80000]; %longer run time to see if equilibrium can be reached 
    end

    for j = 1:3

        if j == 1 %each column is a different density dependnent mortality scenario 
            d2_1 = .05;  % inter < intra
            d1_2 = .05;
        elseif j == 2
            d2_1 = .1; % intra = inter
            d1_2 = .1;
        else
            d2_1 = .15; %inter > intra 
            d1_2 = .15;
        end


        for i = 1:length(propA_vals)
            propA = propA_vals(i);
            envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;

            %first run simulation for partner preference strategy 
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*preference, rtot*preference, rtot*(1-preference), rtot*(1-preference), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            results(1,i) = mean(final_res(1,:));
            results2A(1,i) = median(final_res(1,:));
            results2A(2,i) = quantile(final_res(1,:), .25);
            results2A(3,i) = quantile(final_res(1,:), .75);


            %next run simulation for opposite partner preference strategy 
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*(1-preference), rtot*(1-preference), rtot*preference, rtot*preference, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            results(2,i) = mean(final_res(1,:));

            results2B(1,i) = median(final_res(1,:));
            results2B(2,i) = quantile(final_res(1,:), .25);
            results2B(3,i) = quantile(final_res(1,:), .75);

            %run simulation for reward strategy 50:50
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot/2, rtot/2, rtot/2, rtot/2, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            results(3,i) = mean(final_res(1,:));

            results3(1,i) = median(final_res(1,:));
            results3(2,i) = quantile(final_res(1,:), .25);
            results3(3,i) = quantile(final_res(1,:), .75);

        end

        subplot(2,3,count)
        hold on
        boundedline(propA_vals, results2A(1,:), [results2A(1,:)-results2A(2,:); results2A(3,:)-results2A(1,:)]', 'linewidth', 2, 'color', [196/255 118/255 165/255], 'alpha', 'transparency', 0.5)
        boundedline(propA_vals, results2B(1,:), [results2B(1,:)-results2B(2,:);results2B(3,:)- results2B(1,:)]', 'linewidth', 2, 'color', [128, 180, 232]/255, 'alpha', 'transparency', 0.5)
        boundedline(propA_vals, results3(1,:), [results3(1,:)-results3(2,:);results3(3,:)- results3(1,:)]', 'linewidth', 2, 'color', 'k', 'alpha', 'transparency', 0.5)
        xlabel('Proportion of time in environment A')
        
        count = count+1;

    end
end

%Some labels 
subplot(2,3,1)
title('d_{i,j} < d_{i,i}')

subplot(2,3,2)
title('d_{i,j} = d_{i,i}')
legend({'Bet hedging'; ''; 'F_1 preference'; ''; 'F_2 preference'; ''; ''; ''})

subplot(2,3,3)
title('d_{i,j} > d_{i,i}')
