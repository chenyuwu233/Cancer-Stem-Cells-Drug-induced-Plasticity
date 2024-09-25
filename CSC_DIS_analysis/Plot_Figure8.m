

load('Results\Assumption4_100.mat')
hold on
name = ["\boldmath$p_s$","\boldmath$\alpha_s$","\boldmath$\beta_s$","\boldmath$\nu_{sd}$","\boldmath$p_d$","\boldmath$\alpha_d$"...
"\boldmath$\beta_{d}$","\boldmath$1-b_{d,\beta}$","\boldmath$E_{d,\beta}$","\boldmath$b_{d,\nu}-1$","\boldmath$E_{d,\nu}$"];

% relative_error_a4(:,6) = []; %% temp

sc = boxplot(relative_error_a4,name);
set(sc,'LineWidth',3)
xline([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5])
ax   = gca;
x  = [1,2,3,4,5,6,7,8,9,10,11];
y  = mean(relative_error_a4);
% sa = scatter(x,y,223,"g","x",'LineWidth',3);
set(ax,"TickLabelInterpreter",'latex')
ax.YScale = 'log';
ax.YLim = [1e-4,1e2];
ax.XLim = [0.5,11.5];
yline(0.2,'LineWidth',3)
% yline(0.5,'LineWidth',3)
ylabel('Relative error')
ax.FontSize = 23;
ax.FontWeight = 'bold';