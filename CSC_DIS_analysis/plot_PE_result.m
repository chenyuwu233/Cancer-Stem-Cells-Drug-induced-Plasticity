%% Base case

% load('Results\Base_30.mat')
hold on
% relative_error(:,4) = [];
% relative_difference(:,4) = [];
name = ["\boldmath$\alpha_s$","\boldmath$\beta_s$","\boldmath$\nu_{sd}$",...
"\boldmath$\alpha_{d}$","\boldmath$\beta_{d}$","\boldmath$b_{d,\beta}$","\boldmath$E_{d,\beta}$","\boldmath$b_{d,\nu}$","\boldmath$E_{d,\nu}$"];

sc = boxplot(relative_difference,name);
set(sc,'LineWidth',3)
xline([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5])
x  = [1,2,3,4,5,6,7,8,9];
y  = mean(relative_error);
sa = scatter(x,y,223,"g","x",'LineWidth',3);
ax   = gca;
set(ax,"TickLabelInterpreter",'latex')
% ax.YScale = 'log';
% YL = get(ax,'YTickLabels');
% set(ax,'YTickLabels',YL,'FontWeight','bold')
% ax.YLim = [1e-4,1e2];
ax.YLim = [-5,5];
ax.XLim = [0.5,8.5];
yline(0,'LineWidth',3)
ylabel('Absolute Precentage Error')
ax.FontSize = 23;
ax.FontWeight = 'bold';
title('Base case parameter estimation ')

%% Base case

load('Results\Base_100.mat')
hold on
name = ["\boldmath$\alpha_s$","\boldmath$\beta_s$","\boldmath$\nu_{sd}$","\boldmath$\alpha_d$",...
"\boldmath$\beta_{d}$","\boldmath$1-b_{d,\beta}$","\boldmath$E_{d,\beta}$","\boldmath$b_{d,\nu}-1$","\boldmath$E_{d,\nu}$"];

sc = boxplot(relative_error,name);
set(sc,'LineWidth',3)
xline([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5])
x  = [1,2,3,4,5,6,7,8,9];
y  = mean(relative_error);
% sa = scatter(x,y,223,"g","x",'LineWidth',3);
ax   = gca;
set(ax,"TickLabelInterpreter",'latex')
ax.YScale = 'log';
% YL = get(ax,'YTickLabels');
% set(ax,'YTickLabels',YL,'FontWeight','bold')
ax.YLim = [1e-4,1e2];
ax.XLim = [0.5,9.5];
yline(0.2,'LineWidth',3)
ylabel('Relative Error')
ax.FontSize = 23;
ax.FontWeight = 'bold';
% title('Base case parameter estimation ')


%% Assumption 3

load('Results\Assumption3_100.mat')
hold on
name = ["\boldmath$\alpha_s$","\boldmath$\beta_s$","\boldmath$\nu_{sd}$","\boldmath$b_{s,\beta}$","\boldmath$E_{s,\beta}$","\boldmath$b_{s,\nu}$","\boldmath$E_{s,\nu}$",...
"\boldmath$\alpha_d$","\boldmath$\beta_{d}$","\boldmath$b_{d,\beta}$","\boldmath$E_{d,\beta}$","\boldmath$b_{d,\nu}$","\boldmath$E_{d,\nu}$"];

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



%% Assumption 4

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


%% Hierarchy 
% load('Results\hierarchy.mat')
hold on
target = cap_0_e_hist;
name = ["\boldmath$\alpha_s$","\boldmath$\beta_s$","\boldmath$\nu_{sd}$","\boldmath$\alpha_{d}$",...
"\boldmath$\beta_{d}$","\boldmath$1-b_{d,\beta}$","\boldmath$E_{d,\beta}$","\boldmath$b_{d,\nu}-1$","\boldmath$E_{d,\nu}$"];
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



%%

idx = 10;
MAT = [relative_error(:,idx),relative_error_EP(:,idx)];
name = {'Base','EP'};
plot_vecs(MAT,name,[0,1000],'relative error',' ',[1,0,0;0,0,1],'signed rank')
ax = gca;
ax.YScale = 'log';


