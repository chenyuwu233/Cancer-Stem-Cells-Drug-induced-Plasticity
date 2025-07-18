%%  Setting the Concentration and Time points

init   = 1e4;
s      = 2;
Conc_D = [0,4,8];
NT   = 3;
Time = [0,12,36];
NC   = length(Conc_D);
NR   = 3;
cmd  = 'CSC_DIS';
lam_s_range       = [0,0.1];
lam_d_range       = [0,0.5];



%% Load data (Time)


s      = 2;
Conc_T = [0,4,8];
NT   = 5;
Time_T = [2,6,12,24,48];
NC   = length(Conc_T);
NR   = 3;

Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figures 2A'); % No provided ask if necessary.
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 2B'); % No provided ask if necessary.
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 2B');   % No provided ask if necessary.


DATA_SC_T = [cell2mat(Cell_SC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_SC(3:5,11:16))];
DATA_TC_T = [cell2mat(Cell_TC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_TC(3:5,11:16))];
DATA_DMSO_T = [cell2mat(DMSO_TC(9:11,2:9)),[nan;nan;nan],cell2mat(DMSO_TC(9:11,11:16))];
DATA_DMSO_avg_T = mean(DATA_DMSO_T,1);

DATA_T    = zeros(NT,NC,NR);
DATA_sc_T = zeros(NT,NC,NR);
for i = 1:5
    DATA_T(i,:,:) = DATA_DMSO_avg_T(3*i-2:3*i).*DATA_TC_T(:,3*i-2:3*i)/100;
    DATA_sc_T(i,:,:) = squeeze(DATA_T(i,:,:)).*DATA_SC_T(:,3*i-2:3*i)/100;
end
% for i = 1:5
%     DATA_T(i,:,:) = DATA_TC_T(:,3*i-2:3*i);
%     DATA_sc_T(i,:,:) = DATA_SC_T(:,3*i-2:3*i);
% end

% DATA_T(1,:,:) = Fit_5.a;
% DATA_sc_T(1,:,:) = 0.25/100*Fit_5.a;


%% Plot Time DATA

hold on
% null_DATA = mean(DATA_T(4,1,:));

for i = 1:length(Conc_T)
    DATA_i = squeeze(DATA_T(:,i,:))';
    T_mean = mean(DATA_i,'omitnan');
    T_std  = std(DATA_i,'omitnan');
    errorbar(Time_T,T_mean,T_std,'LineWidth',3)

end

ax = gca;
xlabel('Time','FontSize',23,'FontWeight','bold')
ylabel('Total cell count','FontSize',23,'FontWeight','bold')
ax.FontSize = 23;
ax.FontWeight = 'bold';
xticks(Time_T)
% xticklabels(Time_T)
legend({'0 um','4 um','8 um'})
title("AGS Time experiment")