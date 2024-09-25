%% Load data

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

Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figures 2A');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 2B');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 2B');


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


%% Load data (Conc)

s      = 2;
Conc_D = [0,0.125,0.25,0.5,1,2,4,8,16];
NT   = 2;
Time_D = [0,24];
NC   = length(Conc_D);
NR   = 2;

Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1F');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1E');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 1E');

DATA_SC_D = [cell2mat(Cell_SC(2:10,2:3))];
DATA_TC_D = [cell2mat(Cell_TC(2:10,2:3))];
DMSO_1_D  = cell2mat(DMSO_TC(2,6:10));
DMSO_2_D  = cell2mat(DMSO_TC(2,11:15));


% 
DATA_D    = zeros(NT,NC,NR);
% DATA(1,:,1) = mean(DMSO_1)/4.2931;
% DATA(1,:,2) = mean(DMSO_2)/4.2931; %[0-48]
% DATA(1,:,1) = mean(DMSO_1)/2.4655; % [0-48] rate = 0.0188
% DATA(1,:,2) = mean(DMSO_2)/2.4655;
% DATA(1,:,1) = mean(DMSO_1)/1.9676; % [12-48]
% DATA(1,:,2) = mean(DMSO_2)/1.9676;
DATA_D(1,:,1) = mean(DMSO_1_D)/1.5702; % [24-48]
DATA_D(1,:,2) = mean(DMSO_2_D)/1.5702;

% DATA_D(2,:,1) = DATA_TC_D(:,1)*mean(DMSO_1_D)/100;
% DATA_D(2,:,2) = DATA_TC_D(:,2)*mean(DMSO_2_D)/100;

DATA_D(2,:,1) = DATA_TC_D(:,1);
DATA_D(2,:,2) = DATA_TC_D(:,2);

DATA_sc_D    = zeros(NT,NC,NR);
DATA_sc_D(2,:,1) = DATA_SC_D(:,1);
DATA_sc_D(2,:,2) = DATA_SC_D(:,2);

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


%% Plot Conc DATA

hold on
D_vec = [1,7,8];
null_DATA = mean(DATA_D(1,1,:));
for i = D_vec
    DATA_i = squeeze(DATA_D(:,i,:))';
    T_mean = mean(DATA_i,'omitnan')./null_DATA;
    T_std  = std(DATA_i,'omitnan')./null_DATA;
    errorbar(24+Time_D,T_mean,T_std)
end

ax = gca;
xlabel('Time')
ylabel('Total cell count')
% legend({'0 um','4 um','8 um'})





%% Check the difference

null_DATA_T = mean(DATA_T(5,1,:));
null_DATA_D = mean(DATA_D(2,1,:));

T = squeeze(DATA_T(5,:,:));
D = squeeze(DATA_D(2,:,:));

hold on

TEMP = zeros(9,2);
TEMP(:,1) = mean(D,2);
tp = mean(T,2);
TEMP(1,2) = tp(1);
TEMP(7:8,2) =tp(2:3);

err_low = zeros(9,2);
err_hig = zeros(9,2);
err_low(:,1) = TEMP(:,1)-min(D')';
ll = tp-min(T')';
err_low(1,2) = ll(1);
err_low(7:8,2) =ll(2:3);
err_hig(:,1) = max(D')'-TEMP(:,1);
hh = max(T')'-tp;
err_hig(1,2) = hh(1);
err_hig(7:8,2) =hh(2:3);


vec = [1,2,3,4,5,6,7,8,9];
bar(vec,TEMP)
ax = gca;
ax.XTickLabel = Conc_D;
xlabel('Concentration Levels')
ylabel('Percentage viability')

er = errorbar(vec-0.145,TEMP(:,1),err_low(:,1),err_hig(:,1),'LineWidth',2);   
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er = errorbar(vec+0.145,TEMP(:,2),err_low(:,2),err_hig(:,2),'LineWidth',2);   
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

sig_4 = ttest2(squeeze(DATA_T(5,2,:)),squeeze(DATA_D(2,7,:)))
sig_8 = ttest2(squeeze(DATA_T(5,3,:)),squeeze(DATA_D(2,8,:)))

%% Check the difference

null_DATA_T = mean(DATA_T(5,1,:));
null_DATA_D = mean(DATA_D(2,1,:));

T = squeeze(DATA_sc_T(5,:,:));
D = squeeze(DATA_sc_D(2,:,:));

hold on

TEMP = zeros(9,2);
TEMP(:,1) = mean(D,2);
tp = mean(T,2);
TEMP(1,2) = tp(1);
TEMP(7:8,2) =tp(2:3);

err_low = zeros(9,2);
err_hig = zeros(9,2);
err_low(:,1) = TEMP(:,1)-min(D')';
ll = tp-min(T')';
err_low(1,2) = ll(1);
err_low(7:8,2) =ll(2:3);
err_hig(:,1) = max(D')'-TEMP(:,1);
hh = max(T')'-tp;
err_hig(1,2) = hh(1);
err_hig(7:8,2) =hh(2:3);


vec = [1,2,3,4,5,6,7,8,9];
bar(vec,TEMP)
ax = gca;
ax.XTickLabel = Conc_D;
xlabel('Concentration Levels')
ylabel('CSCs Percentage')

er = errorbar(vec-0.145,TEMP(:,1),err_low(:,1),err_hig(:,1),'LineWidth',2);   
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er = errorbar(vec+0.145,TEMP(:,2),err_low(:,2),err_hig(:,2),'LineWidth',2);   
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

sig_4 = ttest2(squeeze(DATA_T(5,2,:)),squeeze(DATA_D(2,7,:)))
sig_8 = ttest2(squeeze(DATA_T(5,3,:)),squeeze(DATA_D(2,8,:)))

