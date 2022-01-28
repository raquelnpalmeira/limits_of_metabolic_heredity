%% Figure 3 - heatmaps
div = 1; ndays = 150; dt = 1;

%% Parameter sweep of different K_N_C2 and K_N_N
par_number = 15;
lv = linspace(2,3,par_number);
range_par_N = 10.^lv;
%range_par_N = linspace(1,1e05,par_number);
range_par_C2 = fliplr(range_par_N);
%range_par_N = logspace(0,5,par_number);
%range_par_C2 = fliplr(range_par_N);
heatmap_data = zeros(par_number);
tic
parfor i = 1:par_number
    for j = 1:par_number
        [divnperday,mean_N_conc_eq] = intFunctionNucCatC2_parsweep(range_par_C2(i),range_par_N(j),ndays,div);
        heatmap_data_div(i,j) = divnperday;
        heatmap_data_Nconc(i,j) = mean_N_conc_eq;
    end 
    
end
toc

%% Making the heatmaps

figure;
%Heatmap of number of divisions
subplot(1,2,1);
z = surf(round(log10(range_par_N),2),round(log10(range_par_C2),2),heatmap_data_div);
view(0,90);
d = colorbar;
z.EdgeColor = 'none';
xlabel('log K^{N}_{N}');
ylabel('log K_{C_2}^N');
d.Label.String = 'Number of protocell divisions per day at equilibrium';
d.Label.FontSize = 16;
set(gca,'FontSize',14);

%Heatmap of nucleotide concentration
subplot(1,2,2);
s = surf(round(log10(range_par_N),2),round(log10(range_par_C2),2),heatmap_data_Nconc);
view(0,90);
c=colorbar;
s.EdgeColor = 'none';
xlabel('log K^{N}_{N}');
ylabel('log K_{C_2}^N');
c.Label.String = 'Concentration of nucleotides at equilibrium';
c.Label.FontSize = 16;
set(gca,'FontSize',14);

%%%%% If saving is changing the font size - use saveas(1,'Figure2.jpeg');
%%%%% to save the plot
