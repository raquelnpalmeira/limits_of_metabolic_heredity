%% Figure 3a


Km_range = [0,0.33,0.66,0.99];

divn = zeros(4,4);
t = 150;
cost = 0;
parfor i = 1:4
    
   divn(i,1) = intFunctionNucCatFA(Km_range(i),10^(0.9),t,1, cost)
    
    
end

parfor i = 1:4
    
   divn(i,2) = intFunctionNucCatAA(Km_range(i),10^(0.9),t,1, cost)
    
    
end

parfor i = 1:4
    
   divn(i,3) = intFunctionNucCatS(Km_range(i),10^(0.9),t,1, cost)
    
    
end

parfor i = 1:4
    
   divn(i,4) = intFunctionNucCatE(Km_range(i),10^(0.9),t,1, cost)
    
    
end


figure(3);
subplot(2,2,1);
h=plot(divn, 'LineWidth',1.5);
set(h, {'Color'}, {[0.50,0.82,0.95]; [0.73,0.87,0.16]; [1 0.5 0]; [0.54,0.12,0.76]},{'Marker'},{'o';'s';'h';'d'});
yline(divn(1,1), 'r--', 'LineWidth', 1.5);
legend('FA','AA','S','E');
xticks([1,2,3,4]);
xticklabels({'No','Weak','Medium','Strong'});
ylabel('cell divisions per day at equilibrium');
set(gca,'FontSize',14);


%% Figure 3b
divn2 = zeros(4,4);
parfor i = 1:4
    
   divn2(i,1) = intFunctionNucCatFA(Km_range(i),1e3,t,1, cost)
    
    
end

parfor i = 1:4
    
   divn2(i,2) = intFunctionNucCatAA(Km_range(i),1e3,t,1, cost)
    
    
end

parfor i = 1:4
    
   divn2(i,3) = intFunctionNucCatS(Km_range(i),1e3,t,1, cost)
    
    
end

parfor i = 1:4
    
   divn2(i,4) = intFunctionNucCatE(Km_range(i),1e3,t,1, cost)
    
    
end

figure(3)
subplot(2,2,2);
j=plot(divn2, 'LineWidth',1.5);
set(j, {'Color'}, {[0.50,0.82,0.95]; [0.73,0.87,0.16]; [1 0.5 0]; [0.54,0.12,0.76]},{'Marker'},{'o';'s';'h';'d'});
yline(divn2(1,1), 'r--', 'LineWidth', 1.5);
legend('FA','AA','S','E');
xticks([1,2,3,4]);
xticklabels({'No','Weak','Medium','Strong'});
ylabel('cell divisions per day at equilibrium');
set(gca,'FontSize',14);

 %% Figure 3c
divn3 = zeros(4,4);
cost = 1;
parfor i = 1:4
    
   divn3(i,1) = intFunctionNucCatFA(Km_range(i),1e3,t,1, cost)
    
    
end

parfor i = 1:4
    
   divn3(i,2) = intFunctionNucCatAA(Km_range(i),1e3,t,1, cost)
    
    
end

parfor i = 1:4
    
   divn3(i,3) = intFunctionNucCatS(Km_range(i),1e3,t,1, cost)
    
    
end

parfor i = 1:4
    
   divn3(i,4) = intFunctionNucCatE(Km_range(i),1e3,t,1, cost)
    
    
end

figure(3)
subplot(2,2,4);
j=plot(divn3, 'LineWidth',1.5);
set(j, {'Color'}, {[0.50,0.82,0.95]; [0.73,0.87,0.16]; [1 0.5 0]; [0.54,0.12,0.76]},{'Marker'},{'o';'s';'h';'d'});
yline(divn3(1,1), 'r--', 'LineWidth', 1.5);
legend('FA','AA','S','E');
xticks([1,2,3,4]);
xticklabels({'No','Weak','Medium','Strong'});
ylabel('cell divisions per day at equilibrium');
set(gca,'FontSize',14);