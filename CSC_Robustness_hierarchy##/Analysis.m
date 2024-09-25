%% Record the point estimation and its error

ds_GR_hist = [];
da_GR_hist = [];
nu_ds_hist = [];
alpha_d_hist = [];
sp_hist    = [];


for i = 1:30
    name = strcat('Boot_CI_CSC_DIS_pos_',num2str(i),'.mat');
    load(name)
    est_ds_GR = GR50(opt_xx_pe(11),opt_xx_pe(12),max(Conc));
    est_da_GR = GR50(opt_xx_pe(13),opt_xx_pe(14),max(Conc));
    est_sp    = get_stable_p([opt_xx_pe(1)-opt_xx_pe(2),opt_xx_pe(3);opt_xx_pe(10),opt_xx_pe(8)-opt_xx_pe(9)]);
    org_ds_GR = GR50(theta(11),theta(12),max(Conc));
    org_da_GR = GR50(theta(13),theta(14),max(Conc));
    org_sp    = get_stable_p(A);
    ds_GR_hist = [ds_GR_hist,abs(log(est_ds_GR/org_ds_GR))];
    da_GR_hist = [da_GR_hist,abs(log(est_da_GR/org_da_GR))];
    sp_hist   = [sp_hist,abs(log(min(org_sp)/min(est_sp)))];
    nu_ds_hist = [nu_ds_hist, opt_xx_pe(10)];
    alpha_d_hist = [alpha_d_hist, opt_xx_pe(8)];
end

%%
HIST = [alpha_d_hist',nu_ds_hist',sp_hist',ds_GR_hist',da_GR_hist'];


%% estimation Scatter plot 

ax = gca;
ax.FontName = 'Arial';
xticks([1,2])
xticklabels({'\alpha_D','\nu_{DS}'})
ax.FontSize = 25;
ax.FontWeight = 'bold';
yline(0,'LineWidth',2);
ylabel('Point estimation value')
hold on
for i = 1:2
    y = HIST(:,i);
    x = linspace(i-0.25,i+0.25,30);
    s = scatter(ax,x,y,35,'filled','b');
%     s.MarkerFaceAlpha = 0.7;
end

%% Estimation error plot
ax = gca;
ax.FontName = 'Arial';
xticks([3,4,5])
xticklabels({'sp','GR_{Ds}','GR_{Da}'})
ax.FontSize = 25;
ax.FontWeight = 'bold';
yline(0,'LineWidth',2);
ylabel('Absolute log ratio error')
hold on
for i = 3:5
    y = HIST(:,i);
    x = linspace(i-0.25,i+0.25,30);
    s = scatter(ax,x,y,35,'filled','b');
%     s.MarkerFaceAlpha = 0.7;
end


%% Recording all the point estimation

ds_GR_hist = [];
da_GR_hist = [];
nu_ds_hist = [];
alpha_d_hist = [];
sp_hist    = [];
theta_hist = [];
pe_hist    = [];
absolute_error = [];
org_GR = [];

for i = 31:60
    try
        name = strcat('Boot_CI_CSC_DIS_',num2str(i),'.mat');
        load(name)
        est_ds_GR = GR50(opt_xx_pe(11),opt_xx_pe(12),max(Conc));
        est_da_GR = GR50(opt_xx_pe(13),opt_xx_pe(14),max(Conc));
        est_sp    = get_stable_p([opt_xx_pe(1)-opt_xx_pe(2),opt_xx_pe(3);opt_xx_pe(10),opt_xx_pe(8)-opt_xx_pe(9)]);
        org_ds_GR = GR50(theta(11),theta(12),max(Conc));
        org_da_GR = GR50(theta(13),theta(14),max(Conc));
        org_GR    = [org_GR;[org_ds_GR,org_da_GR]];
%         org_sp    = get_stable_p(A);
        ds_GR_hist = [ds_GR_hist,abs(log(est_ds_GR/org_ds_GR))];
        da_GR_hist = [da_GR_hist,abs(log(est_da_GR/org_da_GR))];
%         sp_hist   = [sp_hist,abs(log(min(org_sp)/min(est_sp)))];
        nu_ds_hist = [nu_ds_hist, opt_xx_pe(10)];
        alpha_d_hist = [alpha_d_hist, opt_xx_pe(8)];
        theta_hist = [theta_hist;theta];
        pe_hist    = [pe_hist;opt_xx_pe];
        absolute_error = [absolute_error;abs(opt_xx_pe - theta)];
    catch
    end
end