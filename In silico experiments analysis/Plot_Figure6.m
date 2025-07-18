%%
load('Results\Base_GS_100.mat')

alpha_r_vec = [relative_error(:,1)];
lam_r_vec   = [lam_error(:,1)];
nu_rs_vec   = [relative_error(:,3)];
alpha_s_vec = [relative_error(:,4)];
lam_s_vec   = [lam_error(:,2)];
GR_beta_vec = [GR_hist(:,1)];
b_beta_vec  = [relative_error(:,6)];
E_beta_vec  = [relative_error(:,7)];
GR_nu_vec   = [GR_hist(:,2)];
b_nu_vec    = [relative_error(:,8)];
E_nu_vec    = [relative_error(:,9)];


%% 
load('Results\Base_GS_100_LLN.mat')

alpha_r_vec = [alpha_r_vec,relative_error(:,1)];
lam_r_vec   = [lam_r_vec,lam_error(:,1)];
nu_rs_vec    = [nu_rs_vec,relative_error(:,3)];
alpha_s_vec = [alpha_s_vec,relative_error(:,4)];
lam_s_vec   = [lam_s_vec,lam_error(:,2)];
GR_beta_vec = [GR_beta_vec,GR_hist(:,1)];
b_beta_vec  = [b_beta_vec,relative_error(:,6)];
E_beta_vec  = [E_beta_vec,relative_error(:,7)];
GR_nu_vec   = [GR_nu_vec,GR_hist(:,2)];
b_nu_vec    = [b_nu_vec,relative_error(:,8)];
E_nu_vec    = [E_nu_vec,relative_error(:,9)];


%%  Set colors



%%  Plot comparison


figure
hold on

curr_x = 1;
ax = gca;
All_pos = [];
XTLabels = {};



%%  alpha_r

Posi = [curr_x,curr_x+0.4];
All_pos = [All_pos,curr_x+0.2];
XTLabels{end+1} = "\boldmath$\alpha_r$";
curr_x = curr_x + 1;
bh = boxplot(alpha_r_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);

set(bh, 'LineWidth',3)



sig_cell = {Posi};
sig_vec  = [ranksum(alpha_r_vec(:,1),alpha_r_vec(:,2))];
H = sigstar(sig_cell,sig_vec);



%%  Lam_r

Posi = [curr_x,curr_x+0.4];
All_pos = [All_pos,curr_x+0.2];
XTLabels{end+1} = "\boldmath$\kappa_r$";
curr_x = curr_x + 1;
bh = boxplot(lam_r_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);

set(bh, 'LineWidth',3)



sig_cell = {Posi};
sig_vec  = [ranksum(lam_r_vec(:,1),lam_r_vec(:,2))];
H = sigstar(sig_cell,sig_vec);


%%  nu_r

Posi = [curr_x,curr_x+0.4];
All_pos = [All_pos,curr_x+0.2];
XTLabels{end+1} = "\boldmath$\nu_{rs}$";
curr_x = curr_x + 1;
bh = boxplot(nu_rs_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);

set(bh, 'LineWidth',3)


sig_cell = {Posi};
sig_vec  = [ranksum(nu_rs_vec(:,1),nu_rs_vec(:,2))];
H = sigstar(sig_cell,sig_vec);

%%  alpha_s

Posi = [curr_x,curr_x+0.4];
All_pos = [All_pos,curr_x+0.2];
XTLabels{end+1} = "\boldmath$\alpha_s$";
curr_x = curr_x + 1;
bh = boxplot(alpha_s_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);

set(bh, 'LineWidth',3)


sig_cell = {Posi};
sig_vec  = [ranksum(alpha_s_vec(:,1),alpha_s_vec(:,2))];
H = sigstar(sig_cell,sig_vec);

%%  Lam_s

Posi = [curr_x,curr_x+0.4];
All_pos = [All_pos,curr_x+0.2];
XTLabels{end+1} = "\boldmath$\kappa_s$";
curr_x = curr_x + 1;
bh = boxplot(lam_s_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);

set(bh, 'LineWidth',3)


sig_cell = {Posi};
sig_vec  = [ranksum(lam_s_vec(:,1),lam_s_vec(:,2))];
H = sigstar(sig_cell,sig_vec);


%%  GR_beta

Posi = [curr_x,curr_x+0.4];
All_pos = [All_pos,curr_x+0.2];
XTLabels{end+1} = "\boldmath$GR_{s,\beta}$";
curr_x = curr_x + 1;
bh = boxplot(GR_beta_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);

set(bh, 'LineWidth',3)


sig_cell = {Posi};
sig_vec  = [ranksum(GR_beta_vec(:,1),GR_beta_vec(:,2))];
H = sigstar(sig_cell,sig_vec);


% %%  b_beta
% 
% Posi = [curr_x,curr_x+0.4];
% All_pos = [All_pos,curr_x+0.2];
% XTLabels{end+1} = "\boldmath$b_{s\beta}$";
% curr_x = curr_x + 1;
% boxplot(b_beta_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0])
% 
% 
% sig_cell = {Posi};
% sig_vec  = [ranksum(b_beta_vec(:,1),b_beta_vec(:,2))];
% H = sigstar(sig_cell,sig_vec);
% 
% %%  E_beta
% 
% Posi = [curr_x,curr_x+0.4];
% All_pos = [All_pos,curr_x+0.2];
% XTLabels{end+1} = "\boldmath$E_{s\beta}$";
% curr_x = curr_x + 1;
% boxplot(E_beta_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0])
% 
% 
% sig_cell = {Posi};
% sig_vec  = [ranksum(E_beta_vec(:,1),E_beta_vec(:,2))];
% H = sigstar(sig_cell,sig_vec);

%%  GR_nu

Posi = [curr_x,curr_x+0.4];
All_pos = [All_pos,curr_x+0.2];
XTLabels{end+1} = "\boldmath$GR_{s,\nu}$";
curr_x = curr_x + 1;
bh = boxplot(GR_nu_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);

set(bh, 'LineWidth',3)


sig_cell = {Posi};
sig_vec  = [ranksum(GR_nu_vec(:,1),GR_nu_vec(:,2))];
H = sigstar(sig_cell,sig_vec);


% %%  b_nu
% 
% Posi = [curr_x,curr_x+0.4];
% All_pos = [All_pos,curr_x+0.2];
% XTLabels{end+1} = "\boldmath$b_{s,\nu}$";
% curr_x = curr_x + 1;
% boxplot(b_nu_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0])
% 
% 
% sig_cell = {Posi};
% sig_vec  = [ranksum(b_nu_vec(:,1),b_nu_vec(:,2))];
% H = sigstar(sig_cell,sig_vec);
% 
% 
% %%  E_nu
% 
% Posi = [curr_x,curr_x+0.4];
% All_pos = [All_pos,curr_x+0.2];
% XTLabels{end+1} = "\boldmath$E_{s,\nu}$";
% curr_x = curr_x + 1;
% boxplot(E_nu_vec,'Positions',Posi,'Symbol','o','OutlierSize',3,'Colors',[0,0,0])
% 
% 
% sig_cell = {Posi};
% sig_vec  = [ranksum(E_nu_vec(:,1),E_nu_vec(:,2))];
% H = sigstar(sig_cell,sig_vec);



%% Coloring


GR1_color = [254 129 125];
GR2_color = [129 184 223];
GR1_color = GR1_color./255;
GR2_color = GR2_color./255;
colors     = [GR1_color;GR2_color];

boxObj=findobj(gca,'Tag','Box');
for i=1:length(boxObj)
    if mod(i,2) == 0
        p1 = patch(boxObj(i).XData,boxObj(i).YData,colors(2,:),'FaceAlpha',0.5,'LineWidth',1.1);
    else
        p2 = patch(boxObj(i).XData,boxObj(i).YData,colors(1,:),'FaceAlpha',0.5,'LineWidth',1.1);
    end
end



%% Modify the axis

    
    ax.XLim           = [0.5,curr_x];
    
    yline(0.2,'LineWidth',3)
    ylabel('Relative Error')
    ax.YLim           = [1e-4,1e3];
    ax.YScale         = 'log';
    ax.LineWidth      = 1.1;
    ax.FontSize       = 25;
    ax.FontName       = 'Arial';
    ax.FontWeight     = 'bold';
    ax.XTick          = All_pos;
    ax.XTickLabel     = XTLabels;
    % ax.Title.String   = title;
    ax.Title.FontSize = 27;
    set(ax,"TickLabelInterpreter",'latex')
    % ax.TicklabelInterpreter = 'latex';
    % ax.YLabel.String  = Y_label;
    

    legend([p1,p2],{'CLT Model','LLN Model'},'location','northwest')