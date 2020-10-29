%2 population SIR model 

clc
clear all
close all

UK_all = readmatrix('Compiled_data.xlsx','Sheet','population','Range','A2:B2');
C_UK = readmatrix('Compiled_data.xlsx','Sheet','contacts','Range','A1:B2');

mu_x = 1/1.6;
mu_y = 1/1.6;
N_x = UK_all(1); % x are children, y are adults
N_y = UK_all(2);
R0 = 0.8;
u = 1/1.6;
N_UK = UK_all;
% Calculate total population size:
N_UK_tot = sum(N_UK);
% Calculate sub-population proportions:
N_UK_prop = N_UK/N_UK_tot;
S = transpose(sum(C_UK, 2)); 

T_1 = transpose(bsxfun(@times, C_UK, N_UK_prop));
eig_max_1 = max(eig(T_1));
eig_min_1 = min(eig(T_1));
if abs(eig_min_1) > abs(eig_max_1)
    dom_eig_1 = eig_min_1;
end 
if abs(eig_max_1) > abs(eig_min_1)
    dom_eig_1 = eig_max_1;
end

b = (R0*u)/dom_eig_1; 
    
betaxx = b*C_UK(1,1);
betaxy = b*C_UK(1,2); 
betayx = b*C_UK(2,1);
betayy = b*C_UK(2,2);

S_x = N_x - 1; %we start off with 1 infected in population x, zero in y
I_x = 0;
R_x = 0;

S_y = N_y;
I_y = 1;
R_y = 0;

timings = 0;

tVec = [];
IVec_x = [];
IVec_y = [];

eventNumber = 1;
tVec(eventNumber) = timings;
IVec_x(eventNumber) = I_x;
IVec_y(eventNumber) = I_y;
eventNumber = eventNumber + 1;

%S_x -> I_x either occurs when susceptible from x is infected by infected by
%infectious from x, or when susceptible from x is infected by infectious
%from y, so at rate betaxx*S_x*I_x + betayx*S_x*I_y
%I_x -> R_x occurs when infectious from x recovers, so occurs at rate
% %mu_x*I_x
%S_y -> I_y at rate betayy*S_y*I_y + betaxy*S_y*I_x

while I_x >= 0 && I_y >= 0 %both I_x and I_y must be greater than 0 for the epidemic to occur, and as I_x and I_y are always positive I_x + I_y is a suitable condition. 
    
    totalRate = betaxx*I_x + betaxy*I_x + mu_x*I_x + betayy*I_y + betayx*I_y + mu_y*I_y;
    r = rand();
    timings = timings + (-1/totalRate)*log(r);
    
    %decide what event is (4 ways this can go: infection of x, infection
    % of y, recovery of x, recovery of y)
    
    rtwo = rand();
    
    Bound_1 = (betaxx*I_x + betayx*I_y)/totalRate;
    Bound_2 = (betaxx*S_x*I_x + betayx*I_y + mu_x*I_x)/totalRate;
    Bound_3 = (betaxx*S_x*I_x + betayx*I_y + mu_x*I_x + betayy*I_y + betaxy*I_x)/totalRate;
    
    if rtwo < Bound_1
       eventType = 1; %infection of x
    end
    
    if Bound_1 < rtwo && rtwo < Bound_2
       eventType = 2; %recovery of x
    end
    
    if Bound_2 < rtwo && rtwo < Bound_3
       eventType = 3; %infection of y
    end
    
    if Bound_3 < rtwo
       eventType = 4; %recovery of y
    end
    
    if eventType == 1
       S_x = S_x - 1;
       I_x = I_x + 1;
    end
    
    if eventType == 2
       I_x = I_x - 1;
       R_x = R_x + 1;
    end
    
    if eventType == 3
       S_y = S_y - 1;
       I_y = I_y + 1;
    end
    
    if eventType == 4
       I_y = I_y - 1;
       R_y = R_y + 1;
    end
  
tVec(eventNumber) = timings;
IVec_x(eventNumber) = I_x;
IVec_y(eventNumber) = I_y;
eventNumber = eventNumber + 1;

end

figure 
plot(tVec, IVec_x)
xlabel('Time')
ylabel('Number of infecteds')

hold on
plot(tVec, IVec_y)
xlabel('time')
hold off

legend('population x', 'population y')

