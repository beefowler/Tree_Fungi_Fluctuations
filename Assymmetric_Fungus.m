% If the fungal strains are asymmetrical in their growth efficiencies or mortality rates
% the best strategy continues to be that which sustains equal populations of both fungi,
% although the distribution of rewards that achieves this result is no longer a 50:50 split


% Maybe 3 columns.
% One fungus with low efficiency, moderate efficiency, high effiency
% Could do 3 rows for mortality too
% The other fungus is fixed at our default values

%Ok keep fungus 1 the same. 

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of nonstructural carbon to mychorrhizal carbon pool
s = 0.01; %senesence of nonstructural carbon
lambda = 0.005; %loss of mycorrhizal carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
m1 = 0.005; %fungus 1 mortality
u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B
sN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


u_bar = .5; %mean uptake of Nitrogen by fungus
difference_val = 0.5;
u1_A = u_bar + difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = u_bar -difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = u_bar -difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = u_bar + difference_val;%uptake of Nitrogen by fungus 2 in environment type B

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
tspan = [1 8000]; 
env_period = 365; 

propA_vals = 0:.02:1; 

%initialize results
results = nan(3,length(propA_vals)); 
results2A = nan(3, length(propA_vals)); 
results2B = results2A; 
results3 = results2B; 


%demographic parameters to use  for fungus 2 
e2_vals = [.8*e1 e1 1.25*e1]; 
m2_vals = [.8*m1 m1 1.25*m1];

for c = 1:3 %vary efficiencies of fungus 2
   
    e2 = e2_vals(c); 
 
   for r = 1:3 %vary mortalities of fungus 2 

        m2 = m2_vals(r); 

       for i = 1:length(propA_vals)

            propA = propA_vals(i);
            envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;

           %3 reward strategies each (1)
           %r1 = rtot / (4/5 * (4/5) + 1);
           r1 = rtot / (4/5 + 1);
           r2 = rtot - r1;

           sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, lambda, r1, r1, r2, r2, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
           final_res = deval(sol, tspan(2)-env_period*3:tspan(2));

           results2A(1,i) = median(final_res(1,:));
           results2A(2,i) = quantile(final_res(1,:), .25);
           results2A(3,i) = quantile(final_res(1,:), .75);

           %3 reward strategies each (2)
           swap = r2;
           r2 = r1;
           r1 = swap;

           sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, lambda, r1, r1, r2, r2, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
           final_res = deval(sol, tspan(2)-env_period*3:tspan(2));


           results2B(1,i) = median(final_res(1,:));
           results2B(2,i) = quantile(final_res(1,:), .25);
           results2B(3,i) = quantile(final_res(1,:), .75);

           %3 reward strategies each (3)
           r1 = rtot/2;
           r2 = rtot/2;

           sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, lambda, r1, r1, r2, r2, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
           final_res = deval(sol, tspan(2)-env_period*3:tspan(2));

           results3(1,i) = median(final_res(1,:));
           results3(2,i) = quantile(final_res(1,:), .25);
           results3(3,i) = quantile(final_res(1,:), .75);


       end


           subplot(3,3,(r-1)*3+c)
           hold on
           boundedline(propA_vals, results2A(1,:), [results2A(1,:)-results2A(2,:); results2A(3,:)-results2A(1,:)]', 'linewidth', 2, 'color', [196/255 118/255 165/255], 'alpha', 'transparency', 0.5)
           boundedline(propA_vals, results2B(1,:), [results2B(1,:)-results2B(2,:);results2B(3,:)- results2B(1,:)]', 'linewidth', 2, 'color', [128, 180, 232]/255, 'alpha', 'transparency', 0.5)
           boundedline(propA_vals, results3(1,:), [results3(1,:)-results3(2,:);results3(3,:)- results3(1,:)]', 'linewidth', 2, 'color', 'k', 'alpha', 'transparency', 0.2)
           xlabel('Proportion of time in environment A')
   end

end

%Some labels
subplot(3,3,1)
title('e_{2} < e_{1}')

subplot(3,3,2)
title('e_{2} = e_{1}')

subplot(3,3,3)
title('e_{2} > e_{1}')

%legend({'r_1 > r_2'; ''; 'r_1 = r_2'; '';  'r_2 > r_1'; '';})

subplot(3,3,1)
ylabel('Tree carbon pool')
subplot(3,3,4)
ylabel('Tree carbon pool')
subplot(3,3,7)
ylabel('Tree carbon pool')

