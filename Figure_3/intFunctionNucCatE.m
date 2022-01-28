function divnperday = intFunctionNucCatE(K_N_E,K_C2_N,ndays,div, cost)


tmax = ndays*24*60*60;

%Initial conditions
C2_ini = 0;
E_ini = 0;
FA_ini = 1e8;
AA_1_ini = 1e8;
AA_2_ini = 1e8;
S_ini = 0;
N_ini = 0;

x_ini= [C2_ini E_ini FA_ini AA_1_ini AA_2_ini S_ini N_ini];   
% Simulation
x=zeros(tmax+1,7);
dx_vector = zeros(tmax,7);
x(1,:)=x_ini;
t=0;
tc=0;
for time=1:tmax
    t=t+1;
    dx = partFunctionNucCatE(x(t,:),K_N_E,K_C2_N, cost);
    dx_vector(t,:) = partFunctionNucCatE(x(t,:),K_N_E,K_C2_N, cost);

    
    if find(x(t,:)+round(dx)<0) > 0
        x(t+1,:)=x(t,:)+round(dx);
        x(t+1,find(x(t,:)+round(dx)<0)) = 0;
    else 
        x(t+1,:)=x(t,:)+round(dx);
        
    end
    if div == 1
        if x(t+1,3)>=2e8
        t=t+1;
        disp(['Cell Divided!!',num2str(x(t,3))]);
        x(t+1,1)= (x(t,1))/2;
        x(t+1,2)= (x(t,2))/2;
        x(t+1,3)= (x(t,3))/2;
        x(t+1,4)= (x(t,4))/2;
        x(t+1,5)= (x(t,5))/2;
        x(t+1,6)= (x(t,6))/2;
        x(t+1,7)= (x(t,7))/2;
        
        if  t < 20*60*60*24
            tc = tc+1;
        end
        
        end

    end
        
end
disp (num2str(tmax));
disp (num2str(t));
disp (num2str(tc));
ndays = tmax/(60*60*24);
disp (['n of divisions after 20 days=',num2str((t-tmax)-tc)]);
disp (num2str(((t-tmax)-tc)/((ndays-20))))

%Avogadro's number
AN = 6.02e23;% mol^-1: 

%Size of fatty acid head
fa_head = 2e-17; 
x_eq = x(20*24*3600:end,:);    
cell_vol_eq = zeros(length(x_eq),1);
nuc_conc_eq = zeros(length(x_eq),1);

for i = 1:length(x_eq)
cell_vol_eq(i,1) = ((fa_head*x_eq(i,3))/3) * (sqrt((fa_head*x_eq(i,3))/(4*pi)));
nuc_conc_eq(i,1) = x_eq(i,7)/(AN*cell_vol_eq(i));
end


mean_N_conc_eq = mean(nuc_conc_eq);


divnperday = ((t-tmax)-tc)/((ndays-20));
