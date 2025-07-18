%% Birth

load('Results\hierarchy.mat')

name = ["0","1","2","3"];
hold on
sc = boxplot(birth_hist,name,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);
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

%%

color1 = [129 184 223];
color1 = color1./255;
boxObj=findobj(gca,'Tag','Box');
for i=1:length(boxObj)
    p1 = patch(boxObj(i).XData,boxObj(i).YData,color1,'FaceAlpha',0.5,'LineWidth',1.1);
end