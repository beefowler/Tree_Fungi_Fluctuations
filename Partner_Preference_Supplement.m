% Script to generate Figure S5. 

% Show that partner preference effect is strong even when relative
% difference in reward rate is very small 

clear all
figure

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of nonstructural carbon to mychorrhizal carbon pool
s = 0.01; %senesence of nonstructural carbon
gamma = 0.005; %loss of mycorrhizal carbon pool to environmenta = .04; %allocation of Carbon to Carbon pool
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.005; %fungus 1 mortality
m2 = 0.005; %fungus 2 mortality
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

u_bar = .5; %mean uptake of nitrogen by fungus
difference_val = .5; %severity of environment, high
u1_A = u_bar + difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = u_bar - difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = u_bar - difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = u_bar + difference_val;%uptake of Nitrogen by fungus 2 in environment type B

d1_1 = .1; 
d2_1 = .05; 
d2_2 = .1; 
d1_2 = .05; 

rtot = .8; 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


% Set timespan and environment conditions during timespan 
tspan = [1 30000]; 

%initialize results
propA_vals = 0:.02:1; 
u2_Bvals = 0:.1:1;
results = nan(3,length(propA_vals)); 
results2A = nan(3, length(propA_vals)); 
results2B = results2A; 

env_period = 365; 

figure 
clf 

preference_vals = [0.6 0.52 0.51]; %here are values for each colunn
for d = 1:3;
    preference = preference_vals(d); 

    if d == 3
        tspan = [1 30000]; 
    end

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for partner prefernce strategy 
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, gamma, rtot*preference, rtot*preference, rtot*(1-preference), rtot*(1-preference), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 
    results2A(1,i) = median(final_res(1,:)); 
    results2A(2,i) = quantile(final_res(1,:), .25); 
    results2A(3,i) = quantile(final_res(1,:), .75); 


    if propA == .64;  %choose one simulation to show and plot it 
        subplot(2,3,3+d)
            b = plot(sol.x, sol.y([1 3 4], :));
           % b(1).Color = [10, 156, 0]/255;
            b(1).Color = [89,174,159]/255;
            %b(2).Color = 'k';
            b(2).Color =  [196/255 118/255 165/255];
            b(3).Color = [128, 180, 232]/255;
            %b(5).Color = [214, 71, 90]/255;

            b(1).LineWidth = 2;
            b(2).LineWidth = 2; 
            b(3).LineWidth = 2; 

            b(2).LineStyle = '-.';
            b(3).LineStyle = '--';

            legend({'Tree'; 'Fungus 1'; 'Fungus 2'})
            %title('F_1 preference')

            xticks([0:10*365:30000])
            xticklabels(10*[0:30000/(365)])
            xlabel('Years')
            xlim([70*365 80*365])
    end


    %next run simulation for opposite preference 
    %environments
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, gamma, rtot*(1-preference), rtot*(1-preference), rtot*preference, rtot*preference, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
   final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 

    results2B(1,i) = median(final_res(1,:)); 
    results2B(2,i) = quantile(final_res(1,:), .25); 
    results2B(3,i) = quantile(final_res(1,:), .75); 


    %run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, gamma, rtot/2, rtot/2, rtot/2, rtot/2, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

end

subplot(2,3,d)
plot(propA_vals, results(3,:), 'k', 'linewidth', 2)
hold on 
boundedline(propA_vals, results2A(1,:), [results2A(1,:)-results2A(2,:); results2A(3,:)-results2A(1,:)]', 'linewidth', 2, 'color', [196/255 118/255 165/255], 'alpha', 'transparency', 0.5)
boundedline(propA_vals, results2B(1,:), [results2B(1,:)-results2B(2,:);results2B(3,:)- results2B(1,:)]', 'linewidth', 2, 'color', [128, 180, 232]/255, 'alpha', 'transparency', 0.5)
xlabel('Proportion of time in environment A')
legend({'Bet hedging'; ''; 'F_1 preference'; ''; 'F_2 preference'; ''; ''; ''})

hold on 
scatter(0.64, results2A(1, propA_vals == 0.64), 25, [1 1 1], 'filled', 'MarkerEdgeColor', 'k')

end


%some labels and titles 
subplot(2,3,1) 
title('60:40 preference')
ylabel('Tree carbon pool')
legend off

subplot(2,3,2) 
title('52:48 preference')

subplot(2,3,3)
title('51:49 preference')
legend off

subplot(2,3,4) 
ylabel('Biomass')
legend off

subplot(2,3,6)
legend off

set(findall(gcf,'-property','FontSize'),'FontSize',13)