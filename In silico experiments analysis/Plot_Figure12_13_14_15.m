%% Hierarchy 
load('Results\hierarchy.mat')
hold on
target = cap_3_e_hist; % Plot Figure 12, 13, 14, 15 by setting this value as cap_0_e_hist, cap_1_e_hist, cap_2_e_hist, cap_3_e_hist respectively.
name = ["\boldmath$\alpha_r$","\boldmath$\beta_r$","\boldmath$\nu_{rs}$","\boldmath$\alpha_{s}$",...
"\boldmath$\beta_{s}$","\boldmath$1-b_{s,\beta}$","\boldmath$E_{s,\beta}$","\boldmath$b_{s,\nu}-1$","\boldmath$E_{s,\nu}$"];
sc = boxplot(target,name);
set(sc,'LineWidth',3)
xline([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5])
ax   = gca;
x  = [1,2,3,4,5,6,7,8,9];
y  = mean(target);
% sa = scatter(x,y,223,"g","x",'LineWidth',3);
set(ax,"TickLabelInterpreter",'latex')
% ax.YScale = 'log';
ax.YLim = [-5,5];
ax.XLim = [0.5,9.5];
yline(0,'LineWidth',3)
yline(0.2,'-.','LineWidth',3)
yline(-0.2,'-.','LineWidth',3)
% yline(0.5,'LineWidth',3)
ylabel('Relative difference')
ax.FontSize = 23;
ax.FontWeight = 'bold';
