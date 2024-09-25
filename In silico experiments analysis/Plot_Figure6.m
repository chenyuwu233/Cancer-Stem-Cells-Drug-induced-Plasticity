    
idx = 92;
if idx <= 60
    name = strcat('Results\In silico Base\Boot_CI_CSC_DIS_alt_',num2str(idx),'.mat');
else
    name = strcat('Results\In silico Base\PE_CSC_DIS_',num2str(idx),'.mat');
end
load(name)

est = [opt_xx_pe(1:3),...
        opt_xx_pe(8:9),1-opt_xx_pe(11),opt_xx_pe(12),opt_xx_pe(13)-1,opt_xx_pe(14)];
org = [theta(1:3),...
        theta(8:9),1-theta(11),theta(12),theta(13)-1,theta(14)];
relative_error = abs(est-org)./org;
% e_hist = [e_hist;est-org];
% mean(relative_error);
% ME = [ME,mean(relative_error)];

Theta =  reshape(theta(1:end-1),[],s)';
Theta_est  = reshape(opt_xx_pe(1:end-1),[],s)';

%% Time based plot
t = tiledlayout(1,2);

ax1 = nexttile;


hold on
fig = gcf;
Conc_cell = {'0 \mu m','0.0313 \mu m','0.0625 \mu m','0.125 \mu m','0.25 \mu m','0.375 \mu m','0.5 \mu m','1.25 \mu m','2.5 \mu m','3.75 \mu m','5 \mu m'};
for d = 1:length(Conc)
    DATA_d = squeeze(DATA(:,d,:))';
    DATA_median = median(DATA_d);
    DATA_std    = std(DATA_d)/3;
    errorbar(Time,DATA_median,DATA_std,'LineWidth',3)
end
legend(Conc_cell,'Location','northwest')

ax1.FontSize = 20;
ax1.FontWeight = 'bold';
xlabel('Time')
ylabel('Total cell count')
ntt = strcat('Average RE:', num2str(mean(relative_error)));
title(ntt)

%%
idx = 129;
if idx <= 60
    name = strcat('Results\In silico Base\Boot_CI_CSC_DIS_alt_',num2str(idx),'.mat');
else
    name = strcat('Results\In silico Base\PE_CSC_DIS_',num2str(idx),'.mat');
end
load(name)

est = [opt_xx_pe(1:3),...
        opt_xx_pe(8:9),1-opt_xx_pe(11),opt_xx_pe(12),opt_xx_pe(13)-1,opt_xx_pe(14)];
org = [theta(1:3),...
        theta(8:9),1-theta(11),theta(12),theta(13)-1,theta(14)];
relative_error = abs(est-org)./org;
% e_hist = [e_hist;est-org];
% mean(relative_error);
% ME = [ME,mean(relative_error)];

Theta =  reshape(theta(1:end-1),[],s)';
Theta_est  = reshape(opt_xx_pe(1:end-1),[],s)';

%% Time based plot


ax2 = nexttile;


hold on
fig = gcf;
Conc_cell = {'0 \mu m','0.0313 \mu m','0.0625 \mu m','0.125 \mu m','0.25 \mu m','0.375 \mu m','0.5 \mu m','1.25 \mu m','2.5 \mu m','3.75 \mu m','5 \mu m'};
for d = 1:length(Conc)
    DATA_d = squeeze(DATA(:,d,:))';
    DATA_median = median(DATA_d);
    DATA_std    = std(DATA_d)/3;
    errorbar(Time,DATA_median,DATA_std,'LineWidth',3)
end
legend(Conc_cell,'Location','northwest')

ax2.FontSize = 20;
ax2.FontWeight = 'bold';
ntt = strcat('Average RE:', num2str(mean(relative_error)));
title(ntt)
xlabel('Time')
ylabel('Total cell count')