function dx = partFunctionNucCatFA(x,K_N_FA,K_C2_N, cost)

dx = zeros(1,7);

%Parameters
CO2_sink = 1e-3;
Km = 1e-4;
Km_N = 2;
K_C2 = 1e5;
K_E_N = 1e03;
K_N_N = 0;
Km_N_FA = 6.9*1e-4;

%Concentrations in the sink 
E_sink = 1e-6;
AA_sink = 1e-6;
S_sink = 0;
N_sink = 0; 

%Permeability constants
%p_fa = 1e-10;
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
dt = 1;

%Percentage yields
if K_N_FA <= 0
lambda_C2_FA = 0.376;
else
lambda_C2_FA = 0.376 * (1 + (K_N_FA * ((N/(AN*cell_vol))/((N/(AN*cell_vol))+Km_N_FA))));
end
if lambda_C2_FA > 1 
    lambda_C2_FA = 1;
end
lambda_C2_S = 0.0801*(1-lambda_C2_FA);
lambda_C2_E = 0.016*(1-lambda_C2_FA);
lambda_C2_AA = 0.9039*(1-lambda_C2_FA);
lambda_C2_AA_1 = 0.5*lambda_C2_AA;
lambda_C2_AA_2 = 0.5*lambda_C2_AA;
lambda_S_N = 0.5;
lambda_AA_2_N = 0.5;



delta_C2 =  (AA_1/(AN*cell_vol) * (CO2_sink)/(CO2_sink + Km) * (K_C2 + K_C2 * (1-(K_N_FA * cost)) * K_C2_N * N/(AN*cell_vol)))*dt;%CO2_sink * AA/(AA+Km) *K_C2


delta_E = lambda_C2_E *C2;
delta_FA = ((lambda_C2_FA  *C2)/5);
delta_AA_1 = lambda_C2_AA_1 *C2;
delta_AA_2 = lambda_C2_AA_2 *C2;
delta_S = ((lambda_C2_S *C2)/2);

if S < E && S < AA_2
% DIVIDE THE PRODUCTION OF NUCLEOTIDES INTO CATALYSED AND NOT, SO CAN SUBTRACT FROM E LATER 
delta_N_E = K_E_N * E/(AN*cell_vol) * S/(AN*cell_vol) * AA_2/(AN*cell_vol) * lambda_S_N * S *dt;   
elseif AA_2 < E && AA_2 < S
delta_N_E = K_E_N * E/(AN*cell_vol) * S/(AN*cell_vol) * AA_2/(AN*cell_vol) * lambda_AA_2_N * AA_2 *dt;
else 
delta_N_E = K_E_N * E/(AN*cell_vol) * S/(AN*cell_vol) * AA_2/(AN*cell_vol) * lambda_AA_2_N * E *dt;    
end

if AA_2>S
% DIVIDE THE PRODUCTION OF NUCLEOTIDES INTO CATALYSED AND NOT, SO CAN SUBTRACT FROM E LATER 
delta_N_N =K_N_N * N/(AN*cell_vol) * (S/(AN*cell_vol))/((S/(AN*cell_vol))+Km_N) * (AA_2/(AN*cell_vol))/((AA_2/(AN*cell_vol))+Km_N)* lambda_S_N * S*dt;   %K_N * lambda_S_N * (S/(AN*cell_vol))/((S/(AN*cell_vol))+Km_E) * lambda_AA_2_N * (AA_2/(AN*cell_vol))/((AA_2/(AN*cell_vol))+Km)*E/(AN*cell_vol)*AA_2 *dt;
else
delta_N_N =K_N_N * N/(AN*cell_vol) * (S/(AN*cell_vol))/((S/(AN*cell_vol))+Km_N) * (AA_2/(AN*cell_vol))/((AA_2/(AN*cell_vol))+Km_N)* lambda_AA_2_N * AA_2*dt;   %K_N * lambda_S_N * (S/(AN*cell_vol))/((S/(AN*cell_vol))+Km_E) * lambda_AA_2_N * (AA_2/(AN*cell_vol))/((AA_2/(AN*cell_vol))+Km)*E/(AN*cell_vol)*AA_2 *dt;
end
% 2*CO2 + AA -> C2 / no loss by diffusion since assumed to immediately
% go into E, AA or FA
     
    dx(1) = delta_C2 - delta_E - delta_FA - delta_AA_1 - delta_AA_2 - delta_S;

% C2 -> E (acetyl-phosphate formation from 2-carbon skeleton)
    
    dx(2) = delta_E+(E_sink - E/(AN*cell_vol)) * cell_sa/(cell_vol) * p_e *dt - delta_N_E;%(E_sink - E) * cell_sa/cell_vol * p_e *dt;
  
% 2*C2 -> FA (simple fatty acid formation from 2-carbon skeleton) / no loss
% by diffusion since assumed to immediately embed into the membrane
    
    dx(3) = delta_FA; %+ (FA_sink - FA)/(AN*cell_vol)) * cell_sa/(AN*cell_vol) * p_fa*dt;%(FA_sink - FA) * cell_sa/cell_vol * p_fa*dt;
% C2 -> AA (glycine - simplest amino acid - formation from 2-carbon
% skeleton)

    dx(4) = delta_AA_1+(AA_sink - (AA_1)/(AN*cell_vol))  * cell_sa/(cell_vol) * p_aa*dt;

    dx(5) = delta_AA_2+(AA_sink - (AA_2)/(AN*cell_vol)) * cell_sa/(cell_vol) * p_aa*dt - delta_N_E - delta_N_N;
    
    dx(6) = delta_S + (S_sink - S/(AN*cell_vol)) * cell_sa/(cell_vol) * p_s*dt - delta_N_E - delta_N_N;
    
    dx(7) = delta_N_E +  delta_N_N + (N_sink - N/(AN*cell_vol)) * cell_sa/(cell_vol) * p_n*dt;



 








end
