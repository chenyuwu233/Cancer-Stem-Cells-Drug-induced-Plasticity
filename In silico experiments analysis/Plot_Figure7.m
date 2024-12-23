

load('Results\Assumption3_100.mat')
hold on
name = ["\boldmath$\alpha_r$","\boldmath$\beta_r$","\boldmath$\nu_{rs}$","\boldmath$b_{r,\beta}$","\boldmath$E_{r,\beta}$","\boldmath$b_{r,\nu}$","\boldmath$E_{r,\nu}$",...
"\boldmath$\alpha_s$","\boldmath$\beta_{s}$","\boldmath$b_{s,\beta}$","\boldmath$E_{s,\beta}$","\boldmath$b_{s,\nu}$","\boldmath$E_{s,\nu}$"];

% relative_error_a3(:,8) = []; %% temp

sc = boxplot(relative_error_a3,name);
set(sc,'LineWidth',3)
xline([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.6,10.5,11.5,12.5,13.5])
ax   = gca;
x  = [1,2,3,4,5,6,7,8,9,10,11,12,13];
y  = mean(relative_error_a3);
% sa = scatter(x,y,223,"g","x",'LineWidth',3);
set(ax,"TickLabelInterpreter",'latex')
ax.YScale = 'log';
ax.YLim = [1e-4,1e2];
ax.XLim = [0.5,13.5];
yline(0.2,'LineWidth',3)
% yline(0.5,'LineWidth',3)
ylabel('Relative error')
ax.FontSize = 23;
ax.FontWeight = 'bold';