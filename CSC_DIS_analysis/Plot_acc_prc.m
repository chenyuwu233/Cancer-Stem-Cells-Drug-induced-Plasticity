


load('Results\In silico Base\record30.mat')
%% Initialization

t = tiledlayout(1,3);
% name = ["\textbf{Stable proportion}","\textbf{Cytotoxic $GR_{50}$}","\textbf{Plasticity $GR_{50}$}"];
name = {'Accuracy','Precision'};



% %% Plot Accuracy
% ax1 = nexttile;
% bh1 = boxplot(Accuracy_hist,name);
% ylabel('\textbf{Absolute log ratio}','Interpreter','latex')
% title('\textbf{Accuracy of point estimation}','Interpreter','latex')
% set(bh1,'LineWidth',3)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'TickLabelInterpreter')
% ax1.FontSize = 23;
% ax1.FontWeight = 'bold';
% 
% 
% 
% %% Plot Precision
% 
% ax2 = nexttile;
% bh2 = boxplot(Precision_hist,name);
% ylabel('\textbf{90 percentile CI width}','Interpreter','latex')
% title('\textbf{Precision of estimation}','Interpreter','latex')
% set(bh2,'LineWidth',3)
% set(gca,'TickLabelInterpreter','latex')
% ax2.FontSize = 23;
% ax2.FontWeight = 'bold';

%% Plot Stable proportion
ax1 = nexttile;
bh1 = boxplot([Accuracy_hist(:,1),Precision_hist(:,1)],name);
ylabel('Feasible Range')
title('Stable proportion')
set(bh1,'LineWidth',3)
% set(gca,'TickLabelInterpreter','latex')
ax1.FontSize = 23;
ax1.FontWeight = 'bold';
ax1.YLim = [0,1];

%% Plot Cytotoxic GR50
ax2 = nexttile;
bh2 = boxplot([Accuracy_hist(:,2),Precision_hist(:,2)],name);
ylabel('Feasible Range')
title('Cytotoxic GR_{50}')
set(bh2,'LineWidth',3)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'TickLabelInterpreter')
ax2.FontSize = 23;
ax2.FontWeight = 'bold';
ax2.YLim = [1e-6,1.8301];

%% Plot Plastic GR50
ax3 = nexttile;
bh3 = boxplot([Accuracy_hist(:,3),Precision_hist(:,3)],name);
ylabel('Feasible Range')
title('Plastic GR_{50}')
set(bh3,'LineWidth',3)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'TickLabelInterpreter')
ax3.FontSize = 23;
ax3.FontWeight = 'bold';
ax3.YLim = [1e-6,1.8301];



%% Plot the cv

hold on
name = ["$\alpha_s$","$\beta_s$","$\nu_{sd}$","$\beta_{d}$","$b_{d,\beta}$","$E_{d,\beta}$","$b_{d,\nu}$","$E_{d,\nu}$"];

vec  = ones(1,30);
% x    = [vec,2*vec,3*vec,4*vec,5*vec,6*vec,7*vec,8*vec];
% y    = reshape(cv_hist,[1,240]);
% sc   = plotSpread(relative_error,'xNames',name,'distributionMarkers','o','showMM',3);
sc = boxplot(relative_error,name);
xline([1.5,2.5,3.5,4.5,5.5,6.5,7.5])
% sc1  = scatter(x,y);
% axmarker = sc{3};
ax   = gca;
% ax.XTickLabel = name;
set(ax,"TickLabelInterpreter",'latex')
ax.YScale = 'log';
ax.YLim = [1e-4,1e2];
ax.XLim = [0.5,8.5];
yline(0.2)
ylabel('Relative error')
ax.FontSize = 23;
ax.FontWeight = 'bold';


%% Plot the me

hold on
name = ["$b_{d,\beta}$","$b_{d,\nu}$"];

vec  = ones(1,30);
x    = [vec,2*vec,3*vec,4*vec,5*vec,6*vec,7*vec,8*vec];
y    = reshape(cv_hist,[1,240]);
sc   = plotSpread(maximum_effect,'xNames',name,'distributionMarkers','o');
xline([1.5])
% sc1  = scatter(x,y);
axmarker = sc{3};
ax   = gca;
% ax.XTickLabel = name;
set(ax,"TickLabelInterpreter",'latex')
ax.YScale = 'log';
ax.YLim = [1e-4,1e2];
ax.XLim = [0.5,2.5];
yline(1)
ylabel('Relative error')
ax.FontSize = 23;
ax.FontWeight = 'bold';

