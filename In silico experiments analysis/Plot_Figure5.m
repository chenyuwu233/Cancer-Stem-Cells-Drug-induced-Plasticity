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
    
end



%%  Plot the scatter

%% Theoretical growth rate change

t = tiledlayout(1,2);

ax = nexttile;

scatter(ME,tgrv,'filled')

ax.XScale = 'log';
xlabel('Average parameters estimation relative error')
ylabel('Maximum change in the long-term growth rate')
xline(0.2,'-','RE = 0.2','LineWidth',1.5)
ax.FontWeight = 'bold';


%% Stable proportion change

ax = nexttile;
scatter(ME,spv,'filled')
ax.XScale = 'log';
xlabel('Average parameters estimation relative error')
ylabel('Maximum change in the stable proportion')
xline(0.2,'-','RE = 0.2','LineWidth',1.5)
ax.FontWeight = 'bold';


