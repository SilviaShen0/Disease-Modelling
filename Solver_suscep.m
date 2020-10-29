%% FAST SIMULTANEOUS EQUATION SOLVER with age group-specific susceptibility
% This code reads in a data matrix and numerically solves the system of
% simultaneous equations governing the r_k. It also calculates the 'whole
% population R0' estimate. It then plots a bar graph of the PMOs along with
% a line overlaid at the position of 1-1/R0 for comparison.

%%
% Read in data
C_UK = readmatrix('MUestimates_all_locations_2.xlsx','Sheet','United Kingdom of Great Britain','Range','A1:P16');
UK_all = readmatrix('UK_POP_AGE.xlsx','Range','A2:P2');

%%
% Set parameters
u = 1/1.6; % mu
%R0 = 1.48;

N_UK = 1000*UK_all;
% Calculate total population size:
N_UK_tot = sum(N_UK);
% Calculate sub-population proportions:
N_UK_prop = N_UK/N_UK_tot;
%% 
% create the matrix of transmission rates T 
% define the vector of age group-specific infectiousness (this can change
% with non-pharmaceutical intervention measures)
% Cijs are the values of C_UK
% define the vector of age group-specific susceptibility (K is the vector
% of sigma-values)
K = [0.33, 0.33, 0.33, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.44, 1.44, 1.44];
R0_UK = 3.11;

% Evaluate the values of matrix T, where Tij = Bi x Cij x Kj 

T = bsxfun(@times, C_UK, K);

S = transpose(sum(T, 2));

b = (R0_UK*u)/sum(S.*N_UK_prop);

%% SOLVE!

% Create n symbolic 'r' variables (r1, r2, ... , rn)
syms 'r' [1 16]

% Solve the system of simultaneous equations (defined in function at end of
% script
fun = @(r)myfunc(r,b);
r0=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5];
r=fsolve(fun,r0)

% Compute vector of PMOs
p = 1-r;

% Define vector of age groups
c = categorical({'00-04','05-09', '10-14', '15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79'});

%% PLOT!

figure(2); hold on; box on; set(gca,'Fontsize',14);
plot(1) = bar(c,p,'Facecolor',[0.5 0.5 0.5])
xlabel('Age group of index case','Fontsize',16);
ylabel('Probability of major outbreak (PMO)','Fontsize',16);
plot(2) = yline(1-1/R0_UK, 'color','r','linewidth',2);
PMO = sum(p.*N_UK_prop);
plot(3) = yline(PMO, 'color',[0 0.5 0],'linewidth',2);
legend(plot([2,3]),'1-1/R_0','Weighted average PMO','Fontsize',14);
legend('boxoff')
title('Probability of major outbreak: UK')
ylim([0 0.9])


%% Define simulataneous system

function F = myfunc(r,b)
C = readmatrix('MUestimates_all_locations_2.xlsx','Sheet','United Kingdom of Great Britain','Range','A1:P16');
K = [0.33, 0.33, 0.33, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.44, 1.44, 1.44]; 
T = bsxfun(@times, C, K);
u = 1/1.6;  
R = (b/u)*T;
S = sum(R,2);
    for k = 1:16
        F(k) = -(1+S(k))*r(k) + 1 + r(k)*( R(k,1)*r(1) + R(k,2)*r(2) + R(k,3)*r(3) + R(k,4)*r(4) + R(k,5)*r(5) + R(k,6)*r(6) + R(k,7)*r(7) + R(k,8)*r(8) +  R(k,9)*r(9) + R(k,10)*r(10) + R(k,11)*r(11) + R(k,12)*r(12) + R(k,13)*r(13) + R(k,14)*r(14) + R(k,15)*r(15) + R(k,16)*r(16));
    end
end

