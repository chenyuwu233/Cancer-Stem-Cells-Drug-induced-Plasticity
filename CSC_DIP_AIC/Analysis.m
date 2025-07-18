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
absolute_pe_error = [];
org_GR = [];
relative_error = [];
cv_hist    = [];
maximum_effect = [];


for i = 31:60
%     try
        name = strcat('Result\PE_CSC_DIS_ip_',num2str(i),'.mat');
        load(name)
%         est_ds_GR = GR50(opt_xx_pe(11),opt_xx_pe(12),max(Conc));
%         est_da_GR = GR50(opt_xx_pe(13),opt_xx_pe(14),max(Conc));
%         est_sp    = get_stable_p([opt_xx_pe(1)-opt_xx_pe(2),opt_xx_pe(3);opt_xx_pe(10),opt_xx_pe(8)-opt_xx_pe(9)]);
%         org_ds_GR = GR50(theta(11),theta(12),max(Conc));
%         org_da_GR = GR50(theta(13),theta(14),max(Conc));
%         org_GR    = [org_GR;[org_ds_GR,org_da_GR]];
%         A         = [theta(1)-theta(2),theta(3);theta(10),theta(8)-theta(9)];
%         org_sp    = get_stable_p(A);
%         ds_GR_hist = [ds_GR_hist,abs(log(est_ds_GR/org_ds_GR))];
%         da_GR_hist = [da_GR_hist,abs(log(est_da_GR/org_da_GR))];
%         sp_hist   = [sp_hist,abs(org_sp(1)-est_sp(1))./org_sp(1)];
%         nu_ds_hist = [nu_ds_hist, opt_xx_pe(10)];
%         alpha_d_hist = [alpha_d_hist, opt_xx_pe(8)];
        theta_hist = [theta_hist;theta];
        pe_hist    = [pe_hist;opt_xx_pe];
        est = [opt_xx_pe(1:4),...
                opt_xx_pe(11),1-opt_xx_pe(13),opt_xx_pe(14),opt_xx_pe(15)-1,opt_xx_pe(16)];
        org = [theta(1:4),...
                theta(11),1-theta(13),theta(14),theta(15)-1,theta(16)];
        absolute_pe_error = [absolute_pe_error;[abs(log(est./org))]];
        relative_error = [relative_error;abs(est-org)./org];
%         B_est = [B_parameter(:,1:3),B_parameter(:,9),B_parameter(:,11:14)];
%         B_cv = std(B_est)./mean(B_est);
%         cv_hist = [cv_hist;B_cv];
        me_vec = [abs((1-opt_xx_pe(11))-(1-theta(11)))./(1-theta(11)), abs((opt_xx_pe(13)-1)-(theta(13)-1))./(theta(13)-1)];
        maximum_effect = [maximum_effect;me_vec];
%     catch
%     end
end


%% 

Accuracy_hist = [];
Precision_hist = [];

for i = 31:60
%     try
        name = strcat('Result\Boot_CI_CSC_DIS_',num2str(i),'.mat');
        load(name)
        est_ds_GR = GR50(opt_xx_pe(11),opt_xx_pe(12),max(Conc));
        est_da_GR = GR50(opt_xx_pe(13),opt_xx_pe(14),max(Conc));
        est_sp    = get_stable_p([opt_xx_pe(1)-opt_xx_pe(2),opt_xx_pe(3);opt_xx_pe(10),opt_xx_pe(8)-opt_xx_pe(9)]);
        org_ds_GR = GR50(theta(11),theta(12),max(Conc));
        org_da_GR = GR50(theta(13),theta(14),max(Conc));
        A         = [theta(1)-theta(2),theta(3);theta(10),theta(8)-theta(9)];
        org_sp    = get_stable_p(A);
        est = [est_sp(1),est_ds_GR,est_da_GR];
        org = [org_sp(1),org_ds_GR,org_da_GR];
        Accuracy_hist = [Accuracy_hist;abs(est-org)];
        B_hist = [];
        for j = 1:100
            est_ds_GR = GR50(B_parameter(j,11),B_parameter(j,12),max(Conc));
            est_da_GR = GR50(B_parameter(j,13),B_parameter(j,14),max(Conc));
            est_sp    = get_stable_p([B_parameter(j,1)-B_parameter(j,2),B_parameter(j,3); ...
                                      B_parameter(j,10),B_parameter(j,8)-B_parameter(j,9)]);
            vec = [est_sp(1),est_ds_GR,est_da_GR];
            B_hist    = [B_hist;vec];
        end
        prc_i = prctile(B_hist,[5,95]);
        
        Precision_hist = [Precision_hist;prc_i(2,:)-prc_i(1,:)];

end

%% Save

save('Result\record30.mat',"relative_error","absolute_pe_error","pe_hist","theta_hist")
