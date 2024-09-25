%% Load data
ME = [];
grv = [];
spv = [];
tgrv = [];
e_hist = [];


for idx = 31:130
    % idx = 92;
    
    if idx <= 60
        name = strcat('Results\In silico Base\Boot_CI_CSC_DIS_alt_',num2str(idx),'.mat');
    else
        name = strcat('Results\In silico Base\PE_CSC_DIS_',num2str(idx),'.mat');
    end
    
    load(name)

    % theta(12) 
    % theta(14)
    % Conc
    
    
    est = [opt_xx_pe(1:3),...
            opt_xx_pe(8:9),1-opt_xx_pe(11),opt_xx_pe(12),opt_xx_pe(13)-1,opt_xx_pe(14)];
    org = [theta(1:3),...
            theta(8:9),1-theta(11),theta(12),theta(13)-1,theta(14)];
    relative_error = abs(est-org)./org;
    e_hist = [e_hist;est-org];
    mean(relative_error);
    ME = [ME,mean(relative_error)];
    
    Theta =  reshape(theta(1:end-1),[],s)';
    Theta_est  = reshape(opt_xx_pe(1:end-1),[],s)';
    
    
    if mean(relative_error) > 0.5
        accuracy = 'Figure\(Base) Bad\';
    elseif mean(relative_error) >0.2
        accuracy = 'Figure\(Base) Medium\';
    else
        accuracy = 'Figure\(Base) Good\';
    end
    name1 = strcat(accuracy,'Time plot',num2str(idx),'.jpg');
    name2 = strcat(accuracy,'Moment plot',num2str(idx),'.jpg');
    name3 = strcat(accuracy,'Est plot',num2str(idx),'.jpg');
    name4 = strcat(accuracy,'Est variance dynamic plot',num2str(idx),'.jpg');
    %% Time based plot
    
    % hold on
    % fig = gcf;
    % for d = 1:length(Conc)
    %     DATA_d = squeeze(DATA(:,d,:))';
    %     DATA_median = median(DATA_d);
    %     DATA_std    = std(DATA_d)/3;
    %     errorbar(Time,DATA_median,DATA_std)
    % end
    
    % saveas(fig,name1);
    
    
    % hold off
    % clf
    % %% Dosage based plot
    % 
    % hold on
    % 
    % fig = gcf;
    % ngr = zeros(1,length(Conc));
    % vgr = zeros(1,length(Conc));
    % dgr = zeros(1,length(Conc));
    % dvgr = zeros(1,length(Conc));
    % ti  = Time(2) - Time(1);
    % b    = Theta(:,1)';
    % A    = Drug_A(Theta,0,cmd);
    % p    = get_stable_p(A);
    % init_c = init*p;
    % for i = 1:length(Conc)
    %     Ai = Drug_A(Theta,Conc(i),cmd);
    %     % if i == length(Conc)
    %     %     get_stable_p(Ai)
    %     % end
    %     ngr(i) = max(eig(Ai));
    %     Mi      = get_Mean(Ai,Time);
    %     Sig_i   = zeros(s,s,s);
    %     for j = 1:s
    %         temp = get_sig_j(Ai,b,ti,j);
    %         Sig_i(:,:,j) = temp;
    %     end
    % 
    %     Vari   = get_Cov_mat(init_c,Ai,Sig_i,Time,Mi);
    %     stdi   = sqrt(diag(Vari));
    %     % plot(Time,stdi)
    % 
    %     fi = fit(Time',stdi,'exp1');
    %     vgr(i) = fi.b;
    % 
    %     DATA_i = squeeze(DATA(:,i,:))';
    %     % std(DATA_i)
    %     % pause
    %     fi = fit(Time',mean(DATA_i)','exp1');
    %     dgr(i) = fi.b;
    %     fi = fit(Time',std(DATA_i)','exp1');
    %     % fi.a
    %     dvgr(i) = fi.b;
    %     % plot(fi,Time,std(DATA_i))
    %     % pause
    % end
    % 
    % plot(Conc,ngr,'LineStyle','-.','LineWidth',3)
    % plot(Conc,vgr,'LineStyle','-.','LineWidth',3)
    % plot(Conc,dgr,'LineWidth',3)
    % plot(Conc,dvgr,'LineWidth',3)
    % legend('True Mean','True std','DATA Mean','DATA std')
    % xlabel('Concentration level')
    % ax = gca;
    % ax.YLim = [-0.1 0.2];
    % ax.XScale = 'log';
    % xline(Conc)
    % 
    % saveas(fig,name2);
    % hold off
    % clf
    % clear
    
    %%

    % hold on
    % 
    % fig = gcf;
    tgr = zeros(1,length(Conc));
    gr = zeros(1,length(Conc));
    gr_e = zeros(1,length(Conc));
    dgr = zeros(1,length(Conc));
    ti  = Time(2) - Time(1);
    b    = Theta(:,1)';
    b_e  = Theta_est(:,1)';
    A    = Drug_A(Theta,0,cmd);
    A_e  = Drug_A(Theta_est,0,cmd);
    p    = get_stable_p(A);
    p_e  = get_stable_p(A_e);
    init_c = init*p;
    init_c_e = init*p_e;
    spvi = [];
    for i = 1:length(Conc)
        Ai = Drug_A(Theta,Conc(i),cmd);
        tgr(i) = max(eig(Ai));
        spi = get_stable_p(Ai);
        spvi = [spvi,spi(1)-p(1)];
        Mi = get_Mean(Ai,Time);
        dist = [1000];
        for j = 1:length(Time)-1
            dist = [dist;sum(init_c*Mi(:,:,j))];
        end
        fi = fit(Time',dist,'exp1');
        gr(i) = fi.b;
        Ai_e = Drug_A(Theta_est,Conc(i),cmd);
        Mi_e = get_Mean(Ai_e,Time);
        dist_e = [1000];
        for j = 1:length(Time)-1
            dist_e = [dist_e;sum(init_c_e*Mi_e(:,:,j))];
        end
        fi = fit(Time',dist_e,'exp1');
        gr_e(i) = fi.b;
        
        DATA_i = squeeze(DATA(:,i,:))';
        fi = fit(Time',mean(DATA_i)','exp1');
        % mean(DATA_i)
        dgr(i) = fi.b;
        % pause
    end
    
    tgrv = [tgrv,max(tgr)-min(tgr)];
    spv = [spv,max(spvi)];
    grv = [grv,max(gr)-min(gr)];
    % plot(Conc,gr,'LineStyle','-.','LineWidth',3)
    % plot(Conc,gr_e,'LineStyle','-.','LineWidth',3)
    % plot(Conc,dgr,'LineWidth',3)
    % legend('True Mean','Est Mean','DATA Mean')
    % xlabel('Concentration level')
    % ax = gca;
    % ax.YLim = [-0.1 0.2];
    % ax.XScale = 'log';
    % xline(Conc)
    % saveas(fig,name3);
    % hold off
    % clf
    % clear
    %%
    % hold on
    % 
    % fig = gcf;
    % 
    % vr = zeros(1,length(Conc));
    % vr_e = zeros(1,length(Conc));
    % dvr = zeros(1,length(Conc));
    % ti  = Time(2) - Time(1);
    % b    = Theta(:,1)';
    % b_e  = Theta_est(:,1)';
    % A    = Drug_A(Theta,0,cmd);
    % A_e  = Drug_A(Theta,0,cmd);
    % p    = get_stable_p(A);
    % p_e  = get_stable_p(A_e);
    % init_c = init*p;
    % init_c_e = init*p_e;
    % for i = 1:length(Conc)
    %     Ai = Drug_A(Theta,Conc(i),cmd);
    %     Mi      = get_Mean(Ai,Time);
    %     Sig_i   = zeros(s,s,s);
    %     for j = 1:s
    %         temp = get_sig_j(Ai,b,ti,j);
    %         Sig_i(:,:,j) = temp;
    %     end
    % 
    %     Vari   = get_Cov_mat(init_c,Ai,Sig_i,Time,Mi);
    %     stdi   = sqrt(diag(Vari));
    %     % plot(Time,stdi)
    % 
    %     fi = fit(Time',stdi,'exp1');
    %     fi
    %     vr(i) = fi.b;
    % 
    % 
    %     Ai_e = Drug_A(Theta_est,Conc(i),cmd);
    %     Mi_e      = get_Mean(Ai_e,Time);
    %     Sig_i_e   = zeros(s,s,s);
    %     for j = 1:s
    %         temp = get_sig_j(Ai_e,b_e,ti,j);
    %         Sig_i_e(:,:,j) = temp;
    %     end
    % 
    %     Vari_e   = get_Cov_mat(init_c_e,Ai_e,Sig_i_e,Time,Mi_e);
    %     stdi_e   = sqrt(diag(Vari_e));
    %     % plot(Time,stdi)
    % 
    %     fi_e = fit(Time',stdi_e,'exp1');
    %     vr_e(i) = fi_e.b;
    %     fi_e
    % 
    %     [stdi';stdi_e']
    %     pause
    % 
    %     DATA_i = squeeze(DATA(:,i,:))';
    %     fi_d = fit(Time',std(DATA_i)','exp1');
    %     dvr(i) = fi_d.b;
    % end
    % 
    % plot(Conc,vr,'LineStyle','-.','LineWidth',3)
    % plot(Conc,vr_e,'LineStyle','-.','LineWidth',3)
    % plot(Conc,dvr,'LineWidth',3)
    % legend('True std dynamic','Est std dynamic','DATA std dynamic')
    % xlabel('Concentration level')
    % ax = gca;
    % ax.YLim = [-0.1 0.2];
    % ax.XScale = 'log';
    % xline(Conc)
    % saveas(fig,name4);
    % hold off
    % clf
    % clear
end



%%  Plot the scatter

%% Theoretical growth rate change

t = tiledlayout(1,2);

ax = nexttile;

scatter(ME,tgrv,'filled')
% ax = gca;
ax.XScale = 'log';
xlabel('Average parameters estimation relative error')
ylabel('Maximum change in the long-term growth rate')
xline(0.2,'-','RE = 0.2','LineWidth',1.5)
ax.FontWeight = 'bold';
% ax.FontSize = 23;

%% Stable proportion change

ax = nexttile;
scatter(ME,spv,'filled')
% ax = gca;
ax.XScale = 'log';
xlabel('Average parameters estimation relative error')
ylabel('Maximum change in the stable proportion')
xline(0.2,'-','RE = 0.2','LineWidth',1.5)
ax.FontWeight = 'bold';
% ax.FontSize = 23;

