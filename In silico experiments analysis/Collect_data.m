%% Recording PE (base)

theta_hist = [];
pe_hist    = [];
relative_error = [];
log_ratio_error = [];
dc_lam_error = [];
relative_difference = [];

for i = 31:130
    if i <= 60
        name = strcat('Boot_CI_CSC_DIS_alt_',num2str(i),'.mat');
    else
        name = strcat('PE_CSC_DIS_',num2str(i),'.mat');
    end
    load(name)
    theta_hist = [theta_hist;theta];
    pe_hist    = [pe_hist;opt_xx_pe];
    est = [opt_xx_pe(1:3),...
            opt_xx_pe(8:9),1-opt_xx_pe(11),opt_xx_pe(12),opt_xx_pe(13)-1,opt_xx_pe(14)];
    org = [theta(1:3),...
            theta(8:9),1-theta(11),theta(12),theta(13)-1,theta(14)];
    relative_error = [relative_error;abs(est-org)./org];
    log_ratio_error = [log_ratio_error;est./org];
    dc_lam_error = [dc_lam_error;abs(opt_xx_pe(8)-opt_xx_pe(9))];
    relative_difference = [relative_difference;(est-org)./org];
end


%% Recording PE (base_2)

theta_hist = [];
pe_hist    = [];
relative_error = [];
log_ratio_error = [];
dc_lam_error = [];
relative_difference = [];

for i = 61:130
    if i <= 60
        name = strcat('Boot_CI_CSC_DIS_alt_',num2str(i),'.mat');
    else
        name = strcat('PE_CSC_DIS_',num2str(i),'.mat');
    end
    load(name)
    theta_hist = [theta_hist;theta];
    pe_hist    = [pe_hist;opt_xx_pe];
    est = [opt_xx_pe(1:3),...
            opt_xx_pe(8:9),1-opt_xx_pe(11),opt_xx_pe(12),opt_xx_pe(13)-1,opt_xx_pe(14)];
    org = [theta(1:3),...
            theta(8:9),1-theta(11),theta(12),theta(13)-1,theta(14)];
    relative_error = [relative_error;abs(est-org)./org];
    log_ratio_error = [log_ratio_error;est./org];
    dc_lam_error = [dc_lam_error;abs(opt_xx_pe(8)-opt_xx_pe(9))];
    relative_difference = [relative_difference;(est-org)./org];
end


%% Recording PE (base EP)

theta_hist_EP = [];
pe_hist_EP    = [];
relative_error_EP = [];
log_ratio_error_EP = [];
dc_lam_error_EP = [];

for i = 31:70

    name = strcat('PE_CSC_DIS_EP_',num2str(i),'.mat');

    load(name)
    theta_hist_EP = [theta_hist_EP;theta];
    pe_hist_EP    = [pe_hist_EP;opt_xx_pe];
    est = [opt_xx_pe(1:3),...
            opt_xx_pe(8:9),1-opt_xx_pe(11),opt_xx_pe(12),opt_xx_pe(13)-1,opt_xx_pe(14)];
    org = [theta(1:3),...
            theta(8:9),1-theta(11),theta(12),theta(13)-1,theta(14)];
    relative_error_EP = [relative_error_EP;abs(est-org)./org];
    log_ratio_error_EP = [log_ratio_error_EP;est./org];
    dc_lam_error_EP = [dc_lam_error_EP;abs(opt_xx_pe(8)-opt_xx_pe(9))];
end

%% Recording PE (a3)

theta_hist_a3 = [];
pe_hist_a3    = [];
relative_error_a3 = [];
error_a3 = [];


for i = 31:130
    name = strcat('PE_CSC_DIS_a3_',num2str(i),'.mat');
    load(name)
    theta_hist_a3 = [theta_hist_a3;theta];
    pe_hist_a3    = [pe_hist_a3;opt_xx_pe];
    est = [opt_xx_pe(1:3),opt_xx_pe(4),opt_xx_pe(5),opt_xx_pe(6),opt_xx_pe(7),...
            opt_xx_pe(8:9),opt_xx_pe(11),opt_xx_pe(12),opt_xx_pe(13),opt_xx_pe(14)];
    org = [theta(1:3),theta(4),theta(5),theta(6),theta(7),...
            theta(8:9),theta(11),theta(12),theta(13),theta(14)];
    relative_error_a3 = [relative_error_a3;abs(est-org)./org];
    error_a3 = [error_a3; est - org];
end


%% Recording PE (a4)

theta_hist_a4 = [];
pe_hist_a4    = [];
relative_error_a4 = [];


for i = 31:130
    name = strcat('PE_CSC_DIS_ip_',num2str(i),'.mat');
    load(name)
    theta_hist_a4 = [theta_hist_a4;theta];
    pe_hist_a4    = [pe_hist_a4;opt_xx_pe];
    est = [opt_xx_pe(1:4),...
            opt_xx_pe(9:11),1-opt_xx_pe(13),opt_xx_pe(14),opt_xx_pe(15)-1,opt_xx_pe(16)];
    org = [theta(1:4),...
            theta(9:11),1-theta(13),theta(14),theta(15)-1,theta(16)];
    relative_error_a4 = [relative_error_a4;abs(est-org)./org];
end

%%  sp GR

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

%%  sp GR

Accuracy_alt_hist = [];
Precision_alt_hist = [];

for i = 31:60
%     try
        name = strcat('Result\Boot_CI_CSC_DIS_alt_',num2str(i),'.mat');
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
        Accuracy_alt_hist = [Accuracy_alt_hist;abs(est-org)];
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
        
        Precision_alt_hist = [Precision_alt_hist;prc_i(2,:)-prc_i(1,:)];

end

%% Hierarchy

cap_0_hist = [];
cap_1_hist = [];
cap_2_hist = [];
cap_3_hist = [];
cap_0_e_hist = [];
cap_1_e_hist = [];
cap_2_e_hist = [];
cap_3_e_hist = [];
sp_hist = [];
sp_0_hist = [];
sp_1_hist = [];
sp_2_hist = [];
sp_3_hist = [];
birth_hist = [];

for i =51:92
    try
        name = strcat('PE_CSC_DIS_a12_',num2str(i),'.mat');
        load(name)
        org = [theta(1:3),...
            theta(8:9),1-theta(11),theta(12),theta(13)-1,theta(14)];
        est = [opt_xx_pe_hist(:,1:3),opt_xx_pe_hist(:,8:9),...
            1-opt_xx_pe_hist(:,11),opt_xx_pe_hist(:,12),opt_xx_pe_hist(:,13)-1,opt_xx_pe_hist(:,14)];
        cap_0_hist = [cap_0_hist;abs(est(1,:)-org)./org];
        cap_1_hist = [cap_1_hist;abs(est(2,:)-org)./org];
        cap_2_hist = [cap_2_hist;abs(est(3,:)-org)./org];
        cap_3_hist = [cap_3_hist;abs(est(4,:)-org)./org];
        cap_0_e_hist = [cap_0_e_hist;(est(1,:)-org)./org];
        cap_1_e_hist = [cap_1_e_hist;(est(2,:)-org)./org];
        cap_2_e_hist = [cap_2_e_hist;(est(3,:)-org)./org];
        cap_3_e_hist = [cap_3_e_hist;(est(4,:)-org)./org];
        birth_hist = [birth_hist;opt_xx_pe_hist(:,8)'./theta(8)];

        Ao = [theta(1)-theta(2),theta(3);theta(10),theta(8)-theta(9)];
        A0 = [est(1,1)-est(1,2),est(1,3);0,est(1,8)-est(1,9)];
        A1 = [est(2,1)-est(2,2),est(2,3);0,est(2,8)-est(2,9)];
        A2 = [est(3,1)-est(3,2),est(3,3);0,est(3,8)-est(3,9)];
        A3 = [est(4,1)-est(4,2),est(4,3);0,est(4,8)-est(4,9)];

        spo = get_stable_p(Ao);
        sp0 = get_stable_p(A0);
        sp1 = get_stable_p(A1);
        sp2 = get_stable_p(A2);
        sp3 = get_stable_p(A3);

        sp_hist = [sp_hist;spo(1)];
        sp_0_hist = [sp_0_hist;sp0(1)];
        sp_1_hist = [sp_1_hist;sp1(1)];
        sp_2_hist = [sp_2_hist;sp2(1)];
        sp_3_hist = [sp_3_hist;sp3(1)];

    catch
    end
end


%% Hierarchy

%% Hierarchy

cap_0_hist = [];
cap_1_hist = [];
cap_2_hist = [];
cap_0_e_hist = [];
cap_1_e_hist = [];
cap_2_e_hist = [];
sp_hist = [];
sp_0_hist = [];
sp_1_hist = [];
sp_2_hist = [];
birth_hist = [];

for i =31:60
    % try
        name = strcat('Results\In silico Limited division\Hierarchy_',num2str(i),'.mat');
        load(name)
        org = [theta(1:3),...
            theta(8:9),1-theta(11),theta(12),theta(13)-1,theta(14)];
        est = [opt_xx_pe_asy_hist(:,1:3),opt_xx_pe_asy_hist(:,8:9),...
            1-opt_xx_pe_asy_hist(:,11),opt_xx_pe_asy_hist(:,12),opt_xx_pe_asy_hist(:,13)-1,opt_xx_pe_asy_hist(:,14)];
        cap_0_hist = [cap_0_hist;abs(est(1,:)-org)./org];
        cap_1_hist = [cap_1_hist;abs(est(2,:)-org)./org];
        cap_2_hist = [cap_2_hist;abs(est(3,:)-org)./org];
        cap_0_e_hist = [cap_0_e_hist;(est(1,:)-org)./org];
        cap_1_e_hist = [cap_1_e_hist;(est(2,:)-org)./org];
        cap_2_e_hist = [cap_2_e_hist;(est(3,:)-org)./org];
        birth_hist = [birth_hist;opt_xx_pe_asy_hist(:,8)'./theta(8)];

        Ao = [theta(1)-theta(2),theta(3);theta(10),theta(8)-theta(9)];
        A0 = [est(1,1)-est(1,2),est(1,3);0,est(1,8)-est(1,9)];
        A1 = [est(2,1)-est(2,2),est(2,3);0,est(2,8)-est(2,9)];
        A2 = [est(3,1)-est(3,2),est(3,3);0,est(3,8)-est(3,9)];

        spo = get_stable_p(Ao);
        sp0 = get_stable_p(A0);
        sp1 = get_stable_p(A1);
        sp2 = get_stable_p(A2);

        sp_hist = [sp_hist;spo(1)];
        sp_0_hist = [sp_0_hist;sp0(1)];
        sp_1_hist = [sp_1_hist;sp1(1)];
        sp_2_hist = [sp_2_hist;sp2(1)];

    % catch
    % end
end

%% Save

save('Result\record30.mat',"Precision_hist","Accuracy_hist","cv_hist","relative_error","absolute_pe_error","pe_hist","theta_hist")
