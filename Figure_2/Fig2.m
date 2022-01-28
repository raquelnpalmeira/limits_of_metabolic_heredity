%% Figure 2

vector = zeros(1,11);

K_C2_N = zeros(1,length(0:0.3:3));
power = 0:0.3:3;

for i = 1:length(0:0.3:3)
   
    K_C2_N(i) = 10^(power(i));
   
    
end

for i = 1:11
    
    
    vector(i) = intFunctionNucCat(K_C2_N(i),300,1);
    disp(['Simulation',num2str(i),'done!'])
    
end


fig = figure(4);
semilogx(K_C2_N,vector, 'x', 'Color', 'black', 'MarkerSize',8);
%r = refline(0,0.55);
%set(r,'LineStyle','--','Color','r','LineWidth',1.5);
ylabel('Cell divisions per day at equilibrium');
xlabel('Rate of nucleotide catalysis of carbon fixation (K^{N}_{C_2})');
set(gca,'FontSize',14);