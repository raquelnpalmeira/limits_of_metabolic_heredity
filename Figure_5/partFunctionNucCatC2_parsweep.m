function dx = partFunctionNucCatC2_parsweep(K_C2_N,K_N_N,x)

dx = zeros(1,7);
%THIS COMBINATION OF PARAMETERS WORKS: intFunctionC2(1e-4,1e-2,1e5, 1e4, 1e4, 1e4, 1, 864000,1)

%Parameters
CO2_sink = 1e-3;
Km = 1e-4;
K_C2 = 1e05;
Km_N = 2;
%K_C2_N = 1e03;%1e4;%100;%100;
K_E_N = 1e03;%100;%100;
%K_N_N = 1e03;

%Concentrations in the sink 
E_sink = 1e-6;
AA_sink = 1e-6;
S_sink = 0;
N_sink = 0; 

%Permeability constants
p_aa = 1e-10;
p_e = 1e-10;
p_s = 1e-10;
p_n =1e-10;

%Variables
C2 = x(:,1);
E = x(:,2);
FA = x(:,3);
AA_1 = x(:,4);
AA_2 = x(:,5);
S = x(:,6);
N = x(:,7);

%Avogadro's number
AN = 6.02e23;% mol^-1: 


fa_head = 2e-17;                % dm^2: - fatty acid head group size 0.2 nm^2 from West et al. 

%Cell surface area

cell_sa = fa_head * FA;

%Cell volume
cell_vol = ((fa_head*FA)/3 * sqrt((fa_head*FA)/(4*pi)));% given that each fatty acid has a size = 1 arbitrary unit, 
%then the surface area of the cell would be the number of fatty acids FA
%and therefore the total volume of the cell would be
%FA/3 * sqrt(FA/4*pi)



%Time interval
dt = 1;%second

lambda_C2_S = 0.05;
lambda_C2_E = 0.01;
lambda_C2_FA = 0.376;
lambda_C2_AA = 0.564;
lambda_C2_AA_1 = 0.5*lambda_C2_AA;
lambda_C2_AA_2 = 0.5*lambda_C2_AA;
lambda_AA_2_N = 0.5;
lambda_S_N = 0.5;

%% Dynamic Equations

%C2 synthesis
delta_C2 =  ((AA_1)/(AN*cell_vol)) * ((CO2_sink)/(CO2_sink + Km)) * (K_C2 + K_C2 * K_C2_N * (N/(AN*cell_vol)))*dt;

%E synthesis
delta_E = lambda_C2_E *C2;

%FA synthesis
delta_FA = ((lambda_C2_FA* C2)/5);

%AA_1 synthesis
delta_AA_1 = lambda_C2_AA_1 *C2;

%AA_2 synthesis
delta_AA_2 = lambda_C2_AA_2 *C2;

%S synthesis
delta_S = ((lambda_C2_S *C2)/2);

%N synthesis

if S < E && S < AA_2
% DIVIDE THE PRODUCTION OF NUCLEOTIDES INTO CATALYSED AND NOT, SO CAN SUBTRACT FROM E LATER 
delta_N_E = K_E_N * E/(AN*cell_vol) * S/(AN*cell_vol) * AA_2/(AN*cell_vol) * lambda_S_N * S *dt;   
elseif AA_2 < E && AA_2 < S
delta_N_E = K_E_N * E/(AN*cell_vol) * S/(AN*cell_vol) * AA_2/(AN*cell_vol) * lambda_AA_2_N * AA_2 *dt;
else
delta_N_E = K_E_N * E/(AN*cell_vol) * S/(AN*cell_vol) * AA_2/(AN*cell_vol) * lambda_AA_2_N * E *dt;    
end

if S < E && S < AA_2
% DIVIDE THE PRODUCTION OF NUCLEOTIDES INTO CATALYSED AND NOT, SO CAN SUBTRACT FROM E LATER 
delta_N_N =K_N_N * N/(AN*cell_vol) * (S/(AN*cell_vol))/((S/(AN*cell_vol))+Km_N) * (AA_2/(AN*cell_vol))/((AA_2/(AN*cell_vol))+Km_N)* lambda_S_N * S*dt;   %K_N * lambda_S_N * (S/(AN*cell_vol))/((S/(AN*cell_vol))+Km_E) * lambda_AA_2_N * (AA_2/(AN*cell_vol))/((AA_2/(AN*cell_vol))+Km)*E/(AN*cell_vol)*AA_2 *dt;
elseif AA_2 < E && AA_2 < S
delta_N_N =K_N_N * N/(AN*cell_vol) * (S/(AN*cell_vol))/((S/(AN*cell_vol))+Km_N) * (AA_2/(AN*cell_vol))/((AA_2/(AN*cell_vol))+Km_N)* lambda_AA_2_N * AA_2*dt;   %K_N * lambda_S_N * (S/(AN*cell_vol))/((S/(AN*cell_vol))+Km_E) * lambda_AA_2_N * (AA_2/(AN*cell_vol))/((AA_2/(AN*cell_vol))+Km)*E/(AN*cell_vol)*AA_2 *dt;
else
delta_N_N =K_N_N * N/(AN*cell_vol) * (S/(AN*cell_vol))/((S/(AN*cell_vol))+Km_N) * (AA_2/(AN*cell_vol))/((AA_2/(AN*cell_vol))+Km_N)* lambda_AA_2_N * E *dt;
end

%Total change in number of C2 molecules     
    dx(1) = delta_C2 - delta_E - delta_FA - delta_AA_1 - delta_AA_2 - delta_S;

%Total change in number of E molecules     
    dx(2) = delta_E+(E_sink - E/(AN*cell_vol)) * cell_sa/(cell_vol) * p_e *dt - delta_N_E;
    
%Total change in number of FA molecules 
    dx(3) = delta_FA; % no loss by diffusion since assumed to immediately embed into the membrane

%Total change in number of AA1 molecules    
    dx(4) = delta_AA_1+(AA_sink - (AA_1)/(AN*cell_vol))  * cell_sa/(cell_vol) * p_aa*dt;
       
%Total change in number of AA2 molecules
    dx(5) = delta_AA_2+(AA_sink - (0.1*AA_2)/(AN*cell_vol)) * cell_sa/(cell_vol) * p_aa*dt - delta_N_E - delta_N_N;
    
%Total change in number of S molecules
    dx(6) = delta_S + (S_sink - S/(AN*cell_vol)) * cell_sa/(cell_vol) * p_s*dt - delta_N_E - delta_N_N;
    
%Total change in number of N molecules
    dx(7) = delta_N_E + delta_N_N + (N_sink - N/(AN*cell_vol)) * cell_sa/(cell_vol) * p_n*dt;

end
