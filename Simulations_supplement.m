% Code to show all state variables in a simulation. 
% Generates Supplement Figure S4 

% Constant partner preference strategy first 
% and then merit-based rewards below 
clear all

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
u_bar = .5; %mean uptake of Nitrogen by fungus

%density dependent mortality parameters
d1_1 = .1;
d2_1 = .05;
d2_2 = .1;
d1_2 = .05;

%maximum reward rate
rtot = 0.8;

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

% Set timespan and environment conditions during timespan
tspan = [1 5000];

%environmental regimes
env_period = 365;
propA_vals = 0:.02:1;

%initialize results
results = nan(3,length(propA_vals)); %mean values for each strategy
results2A = nan(3, length(propA_vals)); %median and iqr for preference A
results2B = results2A; %median and iqr for preference B


figure
clf

difference_val = .25
u1_A = u_bar + difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = u_bar -difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = u_bar -difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = u_bar + difference_val;%uptake of Nitrogen by fungus 2 in environment type B

 % Plot one simulation example for each column
propA = .64; %this is the value of alpha in Figure 2, so we'll use that again. 
envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;

%first run simulation for reward strategy 100% fungus 1 in both environments
sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, lambda, rtot, rtot, 0, 0, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
final_res = deval(sol, tspan(2)-env_period*3:tspan(2));


%plot plot plot 
subplot(3,2,1)
yyaxis left
b = plot(sol.x, sol.y([1], :));
b(1).Color = [89,174,159]/255;

b(1).LineWidth = 2;

ax = gca
ax.YColor = [89,174,159]/295;
%ylim([0 10])

%title('Tree Carbon Pools')
ylabel('Nonstructural tree carbon')
yyaxis right
ax2 = gca
ax2.YColor = 'k'
b2 = plot(sol.x, sol.y([2], :));
b2(1).Color = 'k';
b2(1).LineStyle = '--'; 
b2(1).LineWidth = 2;

ylabel('Mycorrhizal carbon')
xticks([0:365:3000])
xticklabels([0:3000/365])
xlabel('Years')
xlim([4*365 8*365])
legend({'Nonstructural carbon'; 'Mycorrhizal carbon'})
title('Reward Fungus 1 only')

subplot(3,2,3)
c = plot(sol.x, sol.y([3 4], :));


c(1).Color =  [196/255 118/255 165/255];
c(2).Color = [128, 180, 232]/255;
c(1).LineWidth = 2;
c(2).LineWidth = 4;
c(1).LineStyle = '--';
c(2).LineStyle = '-.';


legend({'Fungus 1'; 'Fungus 2'})

xticks([0:365:3000])
xticklabels([0:3000/365])
xlabel('Years')
xlim([4*365 8*365])

ylabel('Fungus carbon')
ylim([0 2])

subplot(3,2,5)

b3 = plot(sol.x, sol.y([5], :));
b3(1).Color = [245, 209, 66]/255;

b3(1).LineWidth = 2;

ylim([0 Ntot+1])

hold on

b3 = plot(sol.x, Ntot-sol.y([5], :));
b3(1).Color = 'k';
b3(1).LineStyle = '--'; 

b3(1).LineWidth = 2;
%title('Nitrogen')

ylabel('Nitrogen')
ax = gca
ax.YColor = [245, 209, 66]/295;
ylim([0 10])
legend({'Tree nitrogen'; 'Soil nitrogen'})
xticks([0:365:3000])
xticklabels([0:3000/365])
xlabel('Years')
xlim([4*365 8*365])

set(findall(gcf,'-property','FontSize'),'FontSize',13)

%% Now we want one for a responsive strategy 


% u_bar = .5; %mean uptake of Nitrogen by fungus
% d = 2 %cycle through environmental severity
% difference_vals = [.1 .25 .5];
% difference_val = difference_vals(4-d);
% u1_A = u_bar + difference_val; %uptake of Nitrogen by fungus 1 in environment type A
% u1_B = u_bar -difference_val; %uptake of Nitrogen by fungus 1 in environment type B
% u2_A = u_bar -difference_val; %uptake of Nitrogen by fungus 2 in environment type A
% u2_B = u_bar + difference_val;%uptake of Nitrogen by fungus 2 in environment type B
% 
% propA = .64; % Plot one simulation example for each column
% envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;

%now run simulation for merit-based reward strategy 
sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, lambda, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
final_res = deval(sol, tspan(2)-env_period*3:tspan(2));

%plot in right hand column 

subplot(3,2,2)
yyaxis left
b = plot(sol.x, sol.y([1], :));
b(1).Color = [89,174,159]/255;

b(1).LineWidth = 2;

ax = gca
ax.YColor = [89,174,159]/295;
%ylim([0 10])

%title('Tree Carbon Pools')
ylabel('Nonstructural tree carbon')
yyaxis right
ax2 = gca
ax2.YColor = 'k'
b2 = plot(sol.x, sol.y([2], :));
b2(1).Color = 'k';
b2(1).LineStyle = '--'; 
b2(1).LineWidth = 2;

ylabel('Mycorrhizal carbon')
xticks([0:365:3000])
xticklabels([0:3000/365])
xlabel('Years')
xlim([4*365 8*365])
legend({'Nonstructural carbon'; 'Mycorrhizal carbon'})
title('Merit-based rewards')

subplot(3,2,4)
c = plot(sol.x, sol.y([3 4], :));



c(1).Color =  [196/255 118/255 165/255];
c(2).Color = [128, 180, 232]/255;
c(1).LineWidth = 2;
c(2).LineWidth = 4;
c(1).LineStyle = '--';
c(2).LineStyle = '-.';

%ylim([0 2])

%title('Fungi')
legend({'Fungus 1'; 'Fungus 2'})

xticks([0:365:3000])
xticklabels([0:3000/365])
xlabel('Years')
xlim([4*365 8*365])

ylabel('Fungus carbon')
ylim([0 2])

subplot(3,2,6)

b3 = plot(sol.x, sol.y([5], :));
b3(1).Color = [245, 209, 66]/255;

b3(1).LineWidth = 2;

ylim([0 Ntot+1])

hold on

b3 = plot(sol.x, Ntot-sol.y([5], :));
b3(1).Color = 'k';
b3(1).LineStyle = '--'; 

b3(1).LineWidth = 2;
%title('Nitrogen')

ylabel('Nitrogen')
ax = gca
ax.YColor = [245, 209, 66]/295;
ylim([0 10])
legend({'Tree nitrogen'; 'Soil nitrogen'})
xticks([0:365:3000])
xticklabels([0:3000/365])
xlabel('Years')
xlim([4*365 8*365])

set(findall(gcf,'-property','FontSize'),'FontSize',13)

