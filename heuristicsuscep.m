% A comparison of the stochastic model vs the 'NEW PMO Solver' model WITH SUSCEPTIBILITY:
% we keep the total-population R0 value constant, but vary the beta-values (which
% we calculate using the HEURISTIC METHOD of calculating R0). 

clc
clear all
close all

C_UK = readmatrix('Compiled_data.xlsx','Sheet','contacts','Range','A1:B2');
UK_all = readmatrix('Compiled_data.xlsx','Sheet','population','Range','A2:B2');
N_UK = 1000*UK_all;
n = 100000;
mu_x = 1/1.6;
mu_y = 1/1.6;

R0_max = 8;
R0 = [0:0.5:R0_max];
PMO_stoch_x = zeros(1+2*R0_max,1);
PMO_stoch_y = zeros(1+2*R0_max,1);
PMO_x = zeros(1+2*R0_max,1);
PMO_y = zeros(1+2*R0_max,1);

u = 1/1.6; % mu
N_UK = 1000*UK_all;
% Calculate total population size:
N_UK_tot = sum(N_UK);
% Calculate sub-population proportions:
N_UK_prop = N_UK/N_UK_tot;

K = [0.33, 1]; 
T = bsxfun(@times, C_UK, K);
S = transpose(sum(T, 2)); % S_k values to calculate heuristic value of R0

m_1 = 1;

for j = 0:0.5:8 
    
    %stochastic model for index case in x
counter_x = [];

b_1 = (j*u)/sum(S.*N_UK_prop);

for c = 1:n
    
N_x_1 = N_UK(1);
N_y_1 = N_UK(2);
    
betaxx = b_1*T(1,1);
betaxy = b_1*T(1,2); 
betayx = b_1*T(2,1);
betayy = b_1*T(2,2);

S_x_1 = N_x_1 - 1; %we start off with 1 infected in population x, zero in y
I_x_1 = 1;
R_x_1 = 0;

S_y_1 = N_y_1;
I_y_1 = 0;
R_y_1 = 0;

tVec_1 = [];
IVec_x_1 = [];
IVec_y_1 = [];
    
timings_1 = 0;
eventNumber_1 = 1;
tVec_1(eventNumber_1) = timings_1;
IVec_x_1(eventNumber_1) = I_x_1;
IVec_y_1(eventNumber_1) = I_y_1;
outcome_x = 0;


while (I_x_1 > 0 && I_y_1>=0) || (I_y_1>0 && I_x_1>=0) %both I_x and I_y must be greater than 0 for the epidemic to occur, and as I_x and I_y are always positive I_x + I_y is a suitable condition. 
    
    totalRate = betaxx*I_x_1 + betaxy*I_x_1 + mu_x*I_x_1 + betayy*I_y_1 + betayx*I_y_1 + mu_y*I_y_1;
    r_one = rand();
    timings_1 = timings_1 + (-1/(totalRate))*log(r_one);
    
    r_two = rand();
    
    Bound_1 = (betaxx*I_x_1 + betayx*I_y_1)/totalRate;
    
    Bound_2 = (betaxx*I_x_1 + betayx*I_y_1 +  mu_x*I_x_1)/totalRate;
   
    Bound_3 = (betaxx*I_x_1 + betayx*I_y_1 +  mu_x*I_x_1 + betayy*I_y_1 + betaxy*I_x_1)/totalRate;
            
            
    if r_two < Bound_1
       eventType = 1; %infection of x
    end
    
    if Bound_1 < r_two && r_two < Bound_2
       eventType = 2; %recovery of x
    end
    
    if Bound_2 < r_two && r_two < Bound_3
       eventType = 3; %infection of y
    end
    
    if Bound_3 < r_two
       eventType = 4; %recovery of y
    end
    
    if eventType == 1
       S_x_1 = S_x_1 - 1;
       I_x_1 = I_x_1 + 1;
    end
    
    if eventType == 2
       I_x_1 = I_x_1 - 1;
       R_x_1 = R_x_1 + 1;
    end
    
    if eventType == 3
       S_y_1 = S_y_1 - 1;
       I_y_1 = I_y_1 + 1;
    end
    
    if eventType == 4
       I_y_1 = I_y_1 - 1;
       R_y_1 = R_y_1 + 1;
    end
        
    if I_x_1 + I_y_1 > 40
       outcome_x = 1;
        break         
    end 
        
    tVec_1(eventNumber_1) = timings_1;
    IVec_x_1(eventNumber_1) = I_x_1;
    IVec_y_1(eventNumber_1) = I_y_1;
    
    eventNumber_1 = eventNumber_1 + 1;
    
end

counter_x(c) = outcome_x;

end 

PMO_stoch_x(m_1) = sum(counter_x)/n;

m_1 = m_1+1;

end

m_2 =1;

for j = 0:0.5:8 
    
    % stochastic model for index case in y
counter_y = [];

b_2 = (j*u)/sum(S.*N_UK_prop);

for c = 1:n
    
N_x_1 = N_UK(1);
N_y_1 = N_UK(2);
    
betaxx = b_2*T(1,1);
betaxy = b_2*T(1,2); 
betayx = b_2*T(2,1);
betayy = b_2*T(2,2);

S_x_1 = N_x_1; %we start off with 1 infected in population x, zero in y
I_x_1 = 0;
R_x_1 = 0;

S_y_1 = N_y_1-1;
I_y_1 = 1;
R_y_1 = 0;

tVec_1 = [];
IVec_x_1 = [];
IVec_y_1 = [];
    
timings_1 = 0;
eventNumber_1 = 1;
tVec_1(eventNumber_1) = timings_1;
IVec_x_1(eventNumber_1) = I_x_1;
IVec_y_1(eventNumber_1) = I_y_1;
outcome_y = 0;


while (I_x_1 > 0 && I_y_1>=0) || (I_y_1>0 && I_x_1>=0) %both I_x and I_y must be greater than 0 for the epidemic to occur, and as I_x and I_y are always positive I_x + I_y is a suitable condition. 
    
    totalRate = betaxx*I_x_1 + betaxy*I_x_1 + mu_x*I_x_1 + betayy*I_y_1 + betayx*I_y_1 + mu_y*I_y_1;
    r_one = rand();
    timings_1 = timings_1 + (-1/(totalRate))*log(r_one);
    
    r_two = rand();
    
    Bound_1 = (betaxx*I_x_1 + betayx*I_y_1)/totalRate;
    
    Bound_2 = (betaxx*I_x_1 + betayx*I_y_1 +  mu_x*I_x_1)/totalRate;
   
    Bound_3 = (betaxx*I_x_1 + betayx*I_y_1 +  mu_x*I_x_1 + betayy*I_y_1 + betaxy*I_x_1)/totalRate;
            
            
    if r_two < Bound_1
       eventType = 1; %infection of x
    end
    
    if Bound_1 < r_two && r_two < Bound_2
       eventType = 2; %recovery of x
    end
    
    if Bound_2 < r_two && r_two < Bound_3
       eventType = 3; %infection of y
    end
    
    if Bound_3 < r_two
       eventType = 4; %recovery of y
    end
    
    if eventType == 1
       S_x_1 = S_x_1 - 1;
       I_x_1 = I_x_1 + 1;
    end
    
    if eventType == 2
       I_x_1 = I_x_1 - 1;
       R_x_1 = R_x_1 + 1;
    end
    
    if eventType == 3
       S_y_1 = S_y_1 - 1;
       I_y_1 = I_y_1 + 1;
    end
    
    if eventType == 4
       I_y_1 = I_y_1 - 1;
       R_y_1 = R_y_1 + 1;
    end
        
    if I_x_1 + I_y_1 > 40
       outcome_y = 1;
        break         
    end 
        
    tVec_1(eventNumber_1) = timings_1;
    IVec_x_1(eventNumber_1) = I_x_1;
    IVec_y_1(eventNumber_1) = I_y_1;
    
    eventNumber_1 = eventNumber_1 + 1;
    
end

counter_y(c) = outcome_y;

end 

PMO_stoch_y(m_2) = sum(counter_y)/n;

m_2 = m_2+1;

end

m_3 = 1;

for j = 0:0.5:8
    
b_3 = (j*u)/sum(S.*N_UK_prop);

% Create n symbolic 'r' variables (r1, r2, ... , rn)
syms 'r' [1 2]

% Solve the system of simultaneous equations (defined in function at end of
% script
fun = @(r)myfunc(r,b_3);
r0=[0,0];
r=fsolve(fun,r0);

% Compute vector of PMOs
p = 1-r;

PMO_x(m_3) = p(1);

PMO_y(m_3) = p(2);

m_3 = m_3+1;

end

figure(1); hold on; box on; set(gca,'fontsize',16);

plot(R0, PMO_stoch_x,'color',[0.1 0.8 1],'linewidth',2)
plot(R0, PMO_stoch_y,'color',[1 0 0],'linewidth',2)
plot(R0, PMO_x,'color',[0.1 0.3 0.9],'linestyle','-.','linewidth',2);
plot(R0, PMO_y,'color',[0.5 0 0],'linestyle','-.','linewidth',2);

xlabel('R0 value') 
ylabel('Probability of a major outbreak');
ylim([0 1]);
leg = legend({'stochastic (child)','stochastic (adult)', 'PMO child', 'PMO adult'});
leg.Location = 'southeast'; leg.Box = 'off';


function F = myfunc(r,b_3)
C = readmatrix('Compiled_data.xlsx','Sheet','contacts','Range','A1:B2');
K = [0.33, 1]; 
T = bsxfun(@times, C, K);
u = 1/1.6;  
R = (b_3/u)*T;
S = sum(R,2);
    for k = 1:2
        F(k) = -(1+S(k))*r(k) + 1 + r(k)*( R(k,1)*r(1) + R(k,2)*r(2) );
    end
end

