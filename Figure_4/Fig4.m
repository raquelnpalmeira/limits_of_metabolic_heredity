%% 1


Km_range = [0,0.33,0.66,0.99];

nconc = zeros(4,4);
t = 150;
cost = 0;
parfor i = 1:4
    
   nconc(i,1) = intFunctionNucCatFA_conc(Km_range(i),10^(0.9),t,1, cost)
    
    
end

parfor i = 1:4
    
   nconc(i,2) = intFunctionNucCatAA_conc(Km_range(i),10^(0.9),t,1, cost)
    
    
end

parfor i = 1:4
    
   nconc(i,3) = intFunctionNucCatS_conc(Km_range(i),10^(0.9),t,1, cost)
    
    
end

parfor i = 1:4
    
   nconc(i,4) = intFunctionNucCatE_conc(Km_range(i),10^(0.9),t,1, cost)
    
    
end


figure(3);
subplot(1,2,1);
h=plot(nconc, 'LineWidth',2);
set(h, {'Color'}, {[0.50,0.82,0.95]; [0.73,0.87,0.16]; [1 0.5 0]; [0.54,0.12,0.76]},{'Marker'},{'o';'s';'h';'d'},{'MarkerSize'},{8});
yline(nconc(1,1), 'r--', 'LineWidth', 2);
legend('FA','AA','S','E');
xticks([1,2,3,4]);
xticklabels({'No','Weak','Medium','Strong'});
xlabel('Catalytic strength');
ylabel('cell divisions per day at equilibrium');
title('Weak nucleotide catalysis of CO_2 fixation');
set(gca,'FontSize',14);


%%
nconc2 = zeros(4,4);
parfor i = 1:4
    
   nconc2(i,1) = intFunctionNucCatFA_conc(Km_range(i),1e3,t,1, cost)
    
    
end

parfor i = 1:4
    
   nconc2(i,2) = intFunctionNucCatAA_conc(Km_range(i),1e3,t,1, cost)
    
    
end

parfor i = 1:4
    
   nconc2(i,3) = intFunctionNucCatS_conc(Km_range(i),1e3,t,1, cost)
    
    
end

parfor i = 1:4
    
   nconc2(i,4) = intFunctionNucCatE_conc(Km_range(i),1e3,t,1, cost)
    
    
end

figure(3)
subplot(1,2,2);
j=plot(nconc2, 'LineWidth',2);
set(j, {'Color'}, {[0.50,0.82,0.95]; [0.73,0.87,0.16]; [1 0.5 0]; [0.54,0.12,0.76]},{'Marker'},{'o';'s';'h';'d'},{'MarkerSize'},{8});
yline(nconc2(1,1), 'r--', 'LineWidth', 2);
legend('FA','AA','S','E');
xticks([1,2,3,4]);
xticklabels({'No','Weak','Medium','Strong'});
xlabel('Catalytic strength');
ylabel('Average nucleotide concentration');
title('Strong nucleotide catalysis of CO_2 fixation');
set(gca,'FontSize',14);

%%
% nconc3 = zeros(4,4);
% cost = 1;
% parfor i = 1:4
%     
%    nconc3(i,1) = intFunctionNucCatFA_conc(Km_range(i),1e3,t,1, cost)
%     
%     
% end
% 
% parfor i = 1:4
%     
%    nconc3(i,2) = intFunctionNucCatAA_conc(Km_range(i),1e3,t,1, cost)
%     
%     
% end
% 
% parfor i = 1:4
%     
%    nconc3(i,3) = intFunctionNucCatS_conc(Km_range(i),1e3,t,1, cost)
%     
%     
% end
% 
% parfor i = 1:4
%     
%    nconc3(i,4) = intFunctionNucCatE_conc(Km_range(i),1e3,t,1, cost)
%     
%     
% end
% 
% figure(2)
% subplot(2,2,4);
% j=plot(nconc3, 'LineWidth',2);
% set(j, {'Color'}, {[0.50,0.82,0.95]; [0.73,0.87,0.16]; [1 0.5 0]; [0.54,0.12,0.76]},{'Marker'},{'o';'s';'h';'d'},{'MarkerSize'},{8});
% yline(nconc3(1,1), 'r--', 'LineWidth', 2);
% legend('FA','AA','S','E');
% xticks([1,2,3,4]);
% xticklabels({'No','Weak','Medium','Strong'});
% xlabel('Catalytic strength');
% ylabel('Average nucleotide concentration');
% title('Weak nucleotide catalysis of CO_2 fixation');
% set(gca,'FontSize',14);