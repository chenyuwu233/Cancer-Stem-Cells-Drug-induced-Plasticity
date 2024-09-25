%% Load the data and select the experiment 
%  This script provides code for generating result similar in Figure 3.
%  We have in total 30 experiments indexed from 31-60, while Figure 3 only
%  demonstrates the experiment 53.
%  These index represents different random seed that generate the data and
%  true parameters. 
%  Please use 'id' to looking for correspondent experiment. 
%
%  E.G. id = 53 means experiment 53.

id  = 53;

d_name = strcat('Results/In silico Base/Boot_CI_CSC_DIS_alt_',num2str(id),'.mat');
load(d_name)



%%  Set the color

GR1_color = [254 129 125];
GR2_color = [129 184 223];
GR1_color = GR1_color./255;
GR2_color = GR2_color./255;
Color     = [GR1_color;GR2_color];



%%  Plot the stable proportion (pie)
True_A  = [theta(1)-theta(2),theta(3);theta(10),theta(8)-theta(9)];
Est_A   = [opt_xx_pe(1)-opt_xx_pe(2),opt_xx_pe(3);opt_xx_pe(10),opt_xx_pe(8)-opt_xx_pe(9)];


True_sp = get_stable_p(True_A);
Est_sp  = get_stable_p(Est_A);

t = tiledlayout('flow');
ax1 = nexttile;
h  = pie(True_sp);
set(h(2:2:end),'FontSize',25,'FontWeight','bold');
ax1.Colormap = Color;
title("True stable proportion")
ax2 = nexttile;
h  = pie(Est_sp);
set(h(2:2:end),'FontSize',25,'FontWeight','bold');
ax2.Colormap = Color;
title("Estimated stable proportion")


%% Plot the stable proportion estimation

ax3 = nexttile;


True_sp_line = [];
True_THETA = coeff_reshape(theta(1:end-1),cmd,'Like',2);
for i = 2:length(Conc)
    True_A = Drug_A(True_THETA,Conc(i),cmd);
    True_sp_line = [True_sp_line;get_stable_p(True_A)];
end

PE_sp_line = [];
PE_THETA = coeff_reshape(opt_xx_pe(1:end-1),cmd,'Like',2);
for i = 2:length(Conc)
    PE_A = Drug_A(PE_THETA,Conc(i),cmd);
    PE_sp_line = [PE_sp_line;get_stable_p(PE_A)];
end

True_CSC_sp = True_sp_line(:,1);
True_CNSC_sp = True_sp_line(:,2);

PE_CSC_sp = PE_sp_line(:,1);
PE_CNSC_sp = PE_sp_line(:,2);

Boot_CSC_sp = [];
Boot_CNSC_sp = [];
for i = 1:100
    b_CSC_sp = [];
    b_CNSC_sp = [];
    THETA_i = coeff_reshape(B_parameter(i,1:end-1),cmd,'Like',2);
    for j = 2:length(Conc)
        Est_A = Drug_A(THETA_i,Conc(j),cmd);
        sp = get_stable_p(Est_A);
        b_CSC_sp = [b_CSC_sp, sp(1)];
        b_CNSC_sp = [b_CNSC_sp,sp(2)];
    end
    Boot_CSC_sp = [Boot_CSC_sp;b_CSC_sp];
    Boot_CNSC_sp = [Boot_CNSC_sp;b_CNSC_sp];
end

Prc_CSC_sp = prctile(Boot_CSC_sp,[5,95]);
Prc_CNSC_sp = prctile(Boot_CNSC_sp,[5,95]);


fill_X = [Conc(2:end),flip(Conc(2:end))];
fill_Y = [Prc_CSC_sp(1,:),flip(Prc_CSC_sp(2,:))];

fill(fill_X,fill_Y,GR2_color,'FaceAlpha',0.5,'EdgeColor',GR2_color)

ax3.XScale = 'log';
ax3.XLim = [min(Conc(2:end)),max(Conc(2:end))];
ax3.XTick = Conc(2:end);
ax3.YLim = [0,1];
hold on
plot(Conc(2:end),True_CSC_sp,'-','Color',GR1_color,'LineWidth',3)
% plot(Conc(2:end),PE_CSC_sp,'-','Color',GR2_color,'LineWidth',3)

scatter(Conc(2:end),Boot_CSC_sp,5,'black','filled')
% errorbar(Conc(2:end),PE_CSC_sp,PE_CSC_sp'-Prc_CSC_sp(1,:),Prc_CSC_sp(2,:)-PE_CSC_sp','-o','Color',GR2_color,'LineWidth',3)
ylabel('\textbf{Stable proportion} \boldmath$\pi(d)$','Interpreter','latex')
xlabel('Drug concentration levels d')
title('CSCs stable proportion estimation')

% plot(Conc(2:end),CNSC_sp,'-.o','Color',GR2_color,'LineWidth',3)
% errorbar(Conc(2:end),PE_CNSC_sp,PE_CNSC_sp'-Prc_CNSC_sp(1,:),Prc_CNSC_sp(2,:)-PE_CNSC_sp','-o','Color','b','LineWidth',3)

%% Compute GR

True_GR = [GR50(theta(11),theta(12),max(Conc)),GR50(theta(13),theta(14),max(Conc))];
Est_GR  = [GR50(opt_xx_pe(11),opt_xx_pe(12),max(Conc)),GR50(opt_xx_pe(13),opt_xx_pe(14),max(Conc))];
Est_B_GR = [];
for i = 1:100
    B_GR = [GR50(B_parameter(i,11),B_parameter(i,12),max(Conc)),GR50(B_parameter(i,13),B_parameter(i,14),max(Conc))];
    Est_B_GR = [Est_B_GR;B_GR];
end


%%  Plot GR
ax4 = nexttile;
hold on
y_lim = [0 3];
xlim([0.01 5])
ylim(y_lim)
xticks(Conc)



set(gca,"Xscale","log")
xline(True_GR(1),'-.','LineWidth',3,'Color',GR1_color)
xline(True_GR(2),'-.','LineWidth',3,'Color',GR2_color)


% for i = 1:length(Conc)
%     if Conc(i) < True_GR(1) && Conc(i+1)> True_GR(2)
%         GR1_int = [Conc(i),Conc(i+1)];
%     end
%     if Conc(i) < True_GR(1) && Conc(i+1)>True_GR(2)
%         GR2_int = [Conc(i),Conc(i+1)];
%     end
% end
% 
% GR1_xconf = [GR1_int(1),GR1_int(2),GR1_int(2),GR1_int(1)];
% GR2_xconf = [GR2_int(1),GR2_int(2),GR2_int(2),GR2_int(1)];
% GR1_yconf = [0,0,4,4];
% GR2_yconf = [0,0,4,4];
% 
% fill(GR1_xconf, GR1_yconf, GR1_color,'FaceAlpha',0.3)
% fill(GR2_xconf, GR2_yconf, GR2_color,'FaceAlpha',0.3)

% hl_yconf  = [0.8 0.8 1.2 1.2];
% dyn_yconf = [1.8 1.8 2.2 2.2];
% sto_yconf = [2.8 2.8 3.2 3.2];

% dyn_yconf = [0.8 0.8 1.2 1.2];
% sto_yconf = [1.8 1.8 2.2 2.2];

% fill(GR1_xconf,hl_yconf,GR1_color)
% fill(GR1_xconf,dyn_yconf,GR1_color)
% fill(GR1_xconf,sto_yconf,GR1_color)
% fill(GR2_xconf,hl_yconf,GR2_color)
% fill(GR2_xconf,dyn_yconf,GR2_color)
% fill(GR2_xconf,sto_yconf,GR2_color)

bh1 = boxplot(Est_B_GR,["cytotoxic effect","plasticity effect"],'orientation','horizontal');
% bh1 = boxchart(["cytotoxic effect","plasticity effect"],Est_B_GR,'orientation','horizontal');
set(bh1,'LineWidth',3)



% Color the box

boxObj=findobj(gca,'Tag','Box');
set(boxObj(1),'Color',GR2_color)
set(boxObj(2),'Color',GR1_color)
% 
patch(boxObj(1).XData,boxObj(1).YData,GR2_color,'FaceAlpha',0.5,'LineWidth',1.1,'Edgecolor',GR2_color);
patch(boxObj(2).XData,boxObj(2).YData,GR1_color,'FaceAlpha',0.5,'LineWidth',1.1,'Edgecolor',GR2_color);


% boxplot(GR1',["Phenopop model", "End Points model", "Live Cell Image model"],"orientation","horizontal")
% bh1 = boxplot(GR1',["End-points", "Live cell image"],"orientation","horizontal","Symbol",'*r');

% boxplot(GR2',["Phenopop model", "End Points model", "Live Cell Image model"],"orientation","horizontal")
% bh2 = boxplot(GR2',[ "End-points", "Live cell image"],"orientation","horizontal",'Symbol','*b');
xline(Conc,'LineWidth',1)
xlabel('Drug concentration levels')
title("GR_{50} estimation")
% fontsize(gca,26,"pixels")
xticklabels(Conc)

%%
%  Pie chart
ax4.Layout.Tile = 4;
% ax5.Layout.TileSpan = [1 4];
ax4.Layout.TileSpan = [1 3];

%  Boxplot
% ax4.Layout.Tile = 2;
% ax4.Layout.TileSpan = [1 2];


ax1.FontSize = 18;
ax2.FontSize = 18;% Pie chart
ax3.FontSize = 18;
ax4.FontSize = 18;
ax1.FontWeight = 'bold';
ax2.FontWeight = 'bold';
ax3.FontWeight = 'bold';
ax4.FontWeight = 'bold';
set(bh1,'LineWidth',3)
% set(bh2,'LineWidth',3)



