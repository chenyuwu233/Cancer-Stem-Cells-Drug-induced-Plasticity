%% Birth

load('Results\hierarchy.mat')

name = ["0","1","2","3"];
hold on
sc = boxplot(birth_hist,name);
set(sc,'LineWidth',3)
ax   = gca;
yline(1,'LineWidth',3)
x  = [1,2,3,4];
y  = mean(birth_hist);
% sa = scatter(x,y,223,"g","x",'LineWidth',3);

% ax.YScale = 'log';
ax.YLim = [0,3];
ax.XLim = [0,5];
ylabel('\boldmath$\hat{\alpha_s}/\alpha_s^*$','Interpreter','latex')
xlabel('Maximum division number')

% set(ax,"DefaultAxesTickLabelInterpreter",'latex')
ax.FontSize = 25;
ax.FontWeight = 'bold';