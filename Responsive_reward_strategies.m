%Code to generate Figure 4 in manuscript 


% Responsive strategies with a few differnet values of leakiness. 
% plot tree biomass variance and fungal biomass variance. 
%Mirror image fungi

clear all 
figure 
clf 

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of nonstructural carbon to mychorrhizal carbon pool
s = 0.01; %senesence of nonstructural carbon
gamma = 0.005; %loss of mycorrhizal carbon pool to environment
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

d1_1 = .1; 
d2_1 = .05; 
d2_2 = .1; 
d1_2 = .05; 

% Set timespan and environment conditions during timespan 
tspan = [1 5000]; 
env_period = 365; 
propA_vals = 0:.02:1; 


%initialize results
results = nan(1,length(propA_vals)); %store bet-hedging results
results2 = nan(3, length(propA_vals)); %store responsive strategy results

%want to show fungal biomass as well 
fungi_resultsA = results2; %store fungus 1 responsive strategy results
fungi_resultsB= results2;  %store fungus 3 responsive strategy results
fungi_results = nan(1,length(propA_vals)); %store fungus bet-hedging results (equal for both fungi) 


leakiness_vals = [0 .25 .75]; %three leakiness values to use
rtot = .2;

for col = 1:3 %want each column to be a level of leakiness 
    leakiness = leakiness_vals(col);

for i = 1:length(propA_vals) %vary evenness for  x axis 
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %next run simulation for responsive reward strategy with leak
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, gamma, rtot*(1-leakiness), rtot*leakiness, rtot*leakiness, rtot*(1-leakiness), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));

    results2(1,i) = median(final_res(1,:)); 
    results2(2,i) = quantile(final_res(1,:), .25); 
    results2(3,i) = quantile(final_res(1,:), .75); 

    fungi_resultsA(1,i) = median(final_res(3,:)); 
    fungi_resultsA(2,i) = quantile(final_res(3,:), .25); 
    fungi_resultsA(3,i) = quantile(final_res(3,:), .75); 

    fungi_resultsB(1,i) = median(final_res(4,:)); 
    fungi_resultsB(2,i) = quantile(final_res(4,:), .25); 
    fungi_resultsB(3,i) = quantile(final_res(4,:), .75); 


    %run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, gamma, rtot/2, rtot/2, rtot/2, rtot/2, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(i) = mean(final_res(1,:)); 

    fungi_results(i) = mean(final_res(3,:)); 

end

%Plot tree top row
subplot(2,3,col)
plot(propA_vals, results, 'k', 'linewidth', 2)
hold on 
boundedline(propA_vals, results2(1,:), [results2(1,:)-results2(2,:);results2(3,:)- results2(1,:)]', 'linewidth', 2, 'color', [89,174,159]/255, 'alpha', 'transparency', 0.4)
ylim([0 20])
legend({'Bet-hedging'; ''; 'Switching'})

%Plot fungus bottom row 
subplot(2,3,col+3)
plot(propA_vals, fungi_results, 'k', 'linewidth', 2)
hold on 
boundedline(propA_vals, fungi_resultsA(1,:), [fungi_resultsA(1,:)-fungi_resultsA(2,:); fungi_resultsA(3,:)-fungi_resultsA(1,:)]', 'linewidth', 2, 'color', [196/255 118/255 165/255], 'alpha', 'transparency', 0.5)
boundedline(propA_vals, fungi_resultsB(1,:), [fungi_resultsB(1,:)-fungi_resultsB(2,:);fungi_resultsB(3,:)- fungi_resultsB(1,:)]', '--', 'linewidth', 2, 'color', [128, 180, 232]/255, 'alpha', 'transparency', 0.5)
xlabel('Proportion of time in environment A')
ylim([0 1.8])

end

%some titles
subplot(2,3,1) 
title('No leakiness, "merit-based"')
ylabel('Photosynthetic biomass')

subplot(2,3,2)
title('Slight leakiness (.25)')

subplot(2,3,3)
title('Extreme leakiness (.75)')

subplot(2,3,4) 
ylabel('Fungus carbon pool')

subplot(2,3,6) 
legend({'Bet-hedging: both fungi'; ''; 'Switching: fungus 1'; ''; 'Switching: fungus 2'})


set(findall(gcf,'-property','FontSize'),'FontSize',13)