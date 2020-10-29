% comparing varied infectiousness vs varied susceptibility

clc
clear all
close all

C_UK = readmatrix('Compiled_data.xlsx','Sheet','contacts','Range','A1:B2');
UK_all = readmatrix('Compiled_data.xlsx','Sheet','population','Range','A2:B2');

R0_max = 8;
R0_vec = [0:0.5:R0_max];
PMO_infect_x = zeros(1+2*R0_max,1);
PMO_infect_y = zeros(1+2*R0_max,1);
PMO_suscep_x = zeros(1+2*R0_max,1);
PMO_suscep_y = zeros(1+2*R0_max,1);
PMO_B = zeros(1+2*R0_max,1);
PMO_K = zeros(1+2*R0_max,1);

u = 1/1.6; % mu
N_UK = 1000*UK_all;
% Calculate total population size:
N_UK_tot = sum(N_UK);
% Calculate sub-population proportions:
N_UK_prop = N_UK/N_UK_tot;

M = [0.33, 1]; % either infectiousness, or susceptibility
B = transpose(bsxfun(@times, transpose(C_UK), M)); % transmission rate matrix with varied infectiousness
K = bsxfun(@times, C_UK, M); % transmission rate matrix with varied susceptibility
S_B = transpose(sum(B, 2)); % S_k values to calculate heuristic value of R0
S_K = transpose(sum(K, 2));

m = 1; % varying susceptibility

for j = 0:0.5:R0_max
    
b = (j*u)/sum(S_K.*N_UK_prop);

% Create n symbolic 'r' variables (r1, r2, ... , rn)
syms 'r' [1 2]

% Solve the system of simultaneous equations (defined in function at end of
% script
fun = @(r)myfunc_k(r,b);
r0=[0,0];
r=fsolve(fun,r0);

% Compute vector of PMOs
p = 1-r;

PMO_suscep_x(m) = p(1);

PMO_suscep_y(m) = p(2);

PMO_K(m) = sum(p.*N_UK_prop);

m = m+1;

end

n = 1; % varying infectiousness

for j = 0:0.5:R0_max
    
k = (j*u)/sum(S_B.*N_UK_prop);

% Create n symbolic 'r' variables (r1, r2, ... , rn)
syms 'r' [1 2]

% Solve the system of simultaneous equations (defined in function at end of
% script
fun = @(r)myfunc_b(r,k);
r0=[0,0];
r=fsolve(fun,r0);

% Compute vector of PMOs
p = 1-r;

PMO_infect_x(n) = p(1);

PMO_infect_y(n) = p(2);

PMO_B(n) = sum(p.*N_UK_prop);

n = n+1;

end

figure(1); hold on; box on; set(gca,'fontsize',16);

plot(R0_vec, PMO_infect_x,'color',[0.1 0.8 1],'linewidth',2)
plot(R0_vec, PMO_infect_y,'color',[1 0 0],'linewidth',2)
plot(R0_vec, PMO_suscep_x,'color',[0.1 0.3 0.9],'linestyle','-.','linewidth',2)
plot(R0_vec, PMO_suscep_y,'color',[0.5 0 0],'linestyle','-.','linewidth',2)
plot(R0_vec, PMO_B, '*-k', 'linewidth',2)
plot(R0_vec, PMO_K, 'o-k', 'linewidth',2)


xlabel('R0 value') 
ylabel('Probability of a major outbreak');
ylim([0 1]);
leg = legend({'varying infectiousness (child)','varying infectiousness (adult)', 'varying susceptibility (child)', 'varying susceptibility (adult)', 'Whole population PMO (varying infectiousness)', 'Whole population PMO (varying susceptibility)'});
leg.Location = 'southeast'; leg.Box = 'off';


function F = myfunc_k(r,b)
C = readmatrix('Compiled_data.xlsx','Sheet','contacts','Range','A1:B2');
K = [0.33, 1]; 
T = bsxfun(@times, C, K);
u = 1/1.6;  
R = (b/u)*T;
S = sum(R,2);
    for k = 1:2
        F(k) = -(1+S(k))*r(k) + 1 + r(k)*( R(k,1)*r(1) + R(k,2)*r(2) );
    end
end

function G = myfunc_b(r,k)
C = readmatrix('Compiled_data.xlsx','Sheet','contacts','Range','A1:B2');
B = [0.33, 1]; 
T = transpose(bsxfun(@times, transpose(C), B));
u = 1/1.6;  
R = (k/u)*T;
S = sum(R,2);
    for i = 1:2
        G(i) = -(1+S(i))*r(i) + 1 + r(i)*( R(i,1)*r(1) + R(i,2)*r(2) );
    end
end