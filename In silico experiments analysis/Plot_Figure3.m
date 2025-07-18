%% Base case

load('Results\Base_100.mat')
hold on
name = ["\boldmath$\alpha_r$","\boldmath$\beta_r$","\boldmath$\nu_{rs}$","\boldmath$\alpha_s$",...
"\boldmath$\beta_{s}$","\boldmath$1-b_{s,\beta}$","\boldmath$E_{s,\beta}$","\boldmath$b_{s,\nu}-1$","\boldmath$E_{s,\nu}$"];

sc = boxplot(relative_error,name,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);
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


%% Coloring

color1 = [129 184 223];
color1 = color1./255;
boxObj=findobj(gca,'Tag','Box');
for i=1:length(boxObj)
    p1 = patch(boxObj(i).XData,boxObj(i).YData,color1,'FaceAlpha',0.5,'LineWidth',1.1);
end