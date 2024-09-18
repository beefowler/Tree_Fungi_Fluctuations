% Check if there is an allee effect or if extiction on tree is unstable

clear all 
figure 
 

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of nonstructural carbon to mychorrhizal carbon pool
s = 0.01; %senesence of nonstructural carbon
gamma = 0.005; %loss of mycorrhizal carbon pool to environment
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


%density dependent mortality parameters 
d1_1 = .1; 
d2_1 = .05; 
d2_2 = .1; 
d1_2 = .05; 

% Set timespan and environment conditions during timespan 
tspan = [1 30000]; 

%ALWAYS ENV A for this supplement 
env_period = 365; 
propA = 1; 
envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(4) = 0; %F2
x0(5) = 1; %N

preference = 1;
initial_vals = [0:0.00005:0.001];

rtot_vals = [.1:.05:1]; 

%initialize results
results = nan(length(rtot_vals), length(initial_vals)); 

for j = 1:length(rtot_vals)
    rtot = rtot_vals(j); 

for i = 1:length(initial_vals)
    x0(3) = initial_vals(i); %F1


sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, gamma, rtot*(preference), rtot*(preference), rtot*(1- preference), rtot*(1- preference), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
results(j,i) = mean(final_res(1,:));

end
end

figure
h = pcolor(rtot_vals, initial_vals, (results)');
ylabel('Initial fungus carbon, F(0)')
xlabel('Reward rate, r')
h.EdgeColor = 'none';
colormap('bone')



