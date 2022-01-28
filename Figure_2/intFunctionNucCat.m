function divnperday = intFunctionNucCat(K_C2_N,ndays,div)

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
    dx = partFunctionNucCat(x(t,:),K_C2_N);
    dx_vector(t,:) = partFunctionNucCat(x(t,:),K_C2_N); %create a vector with all dxs

    
    if find(x(t,:)+nearest(dx)<0) > 0
        x(t+1,:)=x(t,:)+nearest(dx);
        x(t+1,find(x(t,:)+nearest(dx)<0)) = 0;
    else 
        x(t+1,:)=x(t,:)+nearest(dx);
        
    end
    if div == 1
        if x(t+1,3)>=2e8
        t=t+1;
        %disp(['Cell Divided!!',num2str(x(t,3))]);
        x(t+1,1)= (x(t,1))/2;
        x(t+1,2)= (x(t,2))/2;
        x(t+1,3)= (x(t,3))/2;
        x(t+1,4)= (x(t,4))/2;
        x(t+1,5)= (x(t,5))/2;
        x(t+1,6)= (x(t,6))/2;
        x(t+1,7)= (x(t,7))/2;
        
        if  t < 10*60*60*24
            tc = tc+1;
        end
        
        end
        
        
        
        
    end
        
end
disp (num2str(tmax));
disp (num2str(t));
disp (num2str(tc));
ndays = tmax/(60*60*24);
disp (['n of divisions after 10 days=',num2str((t-tmax)-tc)]);
disp (num2str(((t-tmax)-tc)/((ndays-10))))

divnperday = ((t-tmax)-tc)/((ndays-10));
