%% Fit cell dynamic

init   = 1e4;
s      = 2;
Conc = [0,4,8];
NT   = 5;
Time = [2,6,12,24,48];
NC   = length(Conc);
NR   = 3;
cmd  = 'CSC_DIS';
lam_s_range       = [0,0.1];
lam_d_range       = [0,0.5];


Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figures 2A');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 2B');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 2B');


DATA_SC = [cell2mat(Cell_SC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_SC(3:5,11:16))];
DATA_TC = [cell2mat(Cell_TC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_TC(3:5,11:16))];
DATA_DMSO = [cell2mat(DMSO_TC(9:11,2:9)),[nan;nan;nan],cell2mat(DMSO_TC(9:11,11:16))];
DATA_DMSO_avg = mean(DATA_DMSO,1);

DATA  = zeros(NT,NC,NR);
DATA_sc = zeros(NT,NC,NR);
for i = 1:5
    k = i;
    DATA(i,:,:) = DATA_DMSO_avg(3*k-2:3*k).*DATA_TC(:,3*k-2:3*k)/100;
    DATA_sc(i,:,:) = squeeze(DATA(i,:,:)).*DATA_SC(:,3*k-2:3*k)/100;
end




E_data = mean(DATA(:,1,:),3,'omitnan');

Fit_5 = fit(Time',E_data,'exp1');


%% Check the early time points growth
Time_e = [2,6,12];
D1_data = mean(DATA(1:3,1,:),3,'omitnan');
% D2_data = mean(DATA(1:3,2,:),3,'omitnan');
% D3_data = mean(DATA(1:3,3,:),3,'omitnan');
Fit_3_D1 = fit(Time_e',D1_data,'exp1');
% Fit_3_D2 = fit(Time_e',D2_data,'exp1');
% Fit_3_D3 = fit(Time_e',D3_data,'exp1');

D1_ng_22 = Fit_3_D1.b;
% D2_ng_22 = Fit_3_D2.b;
% D3_ng_22 = Fit_3_D3.b;
% 
% Drug_effect_22 = [D1_ng_22,D2_ng_22,D3_ng_22] - D1_ng_22;
% 
% Drug_param_22 = Fit_hill_2(Conc,Drug_effect_22);

% init_5T = mean([Fit_3_D1.a,Fit_3_D2.a,Fit_3_D3.a]);
init_5T = Fit_3_D1.a;

ratio = init_5T/mean(DATA(5,1,:));


Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1F');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1E');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 1E');

DATA_SC = [cell2mat(Cell_SC(2:10,2:3))];
DATA_TC = [cell2mat(Cell_TC(2:10,2:3))];
DMSO_1  = cell2mat(DMSO_TC(2,6:10));
DMSO_2  = cell2mat(DMSO_TC(2,11:15));


% init_2T = mean([DMSO_1,DMSO_2])*ratio;
init_2T = [mean(DMSO_1),mean(DMSO_2)]*ratio;

T12_2T = init_2T*exp(D1_ng_22*12);

% %% Fit 24-48 Drug effect
% 
% Time_1 = [24,48];
% D1_data = mean(DATA(4:5,1,:),3,'omitnan');
% D2_data = mean(DATA(4:5,2,:),3,'omitnan');
% D3_data = mean(DATA(4:5,3,:),3,'omitnan');
% % D1_data = squeeze(DATA(4:5,1,:));
% % D2_data = squeeze(DATA(4:5,2,:));
% % D3_data = squeeze(DATA(4:5,3,:));
% Fit_2_D1 = fit(Time_1',D1_data,'exp1');
% Fit_2_D2 = fit(Time_1',D2_data,'exp1');
% Fit_2_D3 = fit(Time_1',D3_data,'exp1');
% 
% 
% D1_ng_24 = Fit_2_D1.b;
% D2_ng_24 = Fit_2_D2.b;
% D3_ng_24 = Fit_2_D3.b;
% 
% Drug_effect_24 = [D1_ng_24,D2_ng_24,D3_ng_24] - D1_ng_24;
% 
% Drug_param_24 = Fit_hill_2(Conc,Drug_effect_24);

% Drug_param_24 = Fit_hill_3(Conc,Drug_effect_24);


%% Fit 12-24 Drug effect

Time_2 = [12,24];
D1_data = mean(DATA(3:4,1,:),3,'omitnan');
D2_data = mean(DATA(3:4,2,:),3,'omitnan');
D3_data = mean(DATA(3:4,3,:),3,'omitnan');
Fit_2_D1 = fit(Time_2',D1_data,'exp1');
Fit_2_D2 = fit(Time_2',D2_data,'exp1');
Fit_2_D3 = fit(Time_2',D3_data,'exp1');


D1_ng_12 = Fit_2_D1.b;
D2_ng_12 = Fit_2_D2.b;
D3_ng_12 = Fit_2_D3.b;

Drug_effect_12 = [D1_ng_12,D2_ng_12,D3_ng_12] - D1_ng_12;

Drug_param_12 = Fit_hill_2(Conc,Drug_effect_12)


%% Fit 12-48 Drug effect

Time_3 = [12,24,48];
D1_data = mean(DATA(3:5,1,:),3,'omitnan');
D2_data = mean(DATA(3:5,2,:),3,'omitnan');
D3_data = mean(DATA(3:5,3,:),3,'omitnan');
Fit_3_D1 = fit(Time_3',D1_data,'exp1');
Fit_3_D2 = fit(Time_3',D2_data,'exp1');
Fit_3_D3 = fit(Time_3',D3_data,'exp1');


D1_ng_18 = Fit_3_D1.b;
D2_ng_18 = Fit_3_D2.b;
D3_ng_18 = Fit_3_D3.b;

Drug_effect_18 = [D1_ng_18,D2_ng_18,D3_ng_18] - D1_ng_18;

Drug_param_18 = Fit_hill_2(Conc,Drug_effect_18);

rate_12 = exp(D1_ng_18*36);



% %% Estimate the spectrum of drug affect growth rate between 24-48 hour
% 
% 
% Conc = [0,0.125,0.25,0.5,1,2,4,8,16];
% 
% Drug_effect = [];
% for i = Conc
%     Drug_effect = [Drug_effect,log(get_Hill2(i,Drug_param_24(1),Drug_param_24(2)))];
% end
% 
% ng = D1_ng_24 + Drug_effect;
% rate = exp(ng.*24);


%% Estimate the spectrum of drug affect growth rate between 12-24 hour

Conc2 = [0,0.125,0.25,0.5,1,2,4,8,16];
Drug_effect = [];
for i = Conc2
    Drug_effect = [Drug_effect,log(get_Hill2(i,Drug_param_12(1),Drug_param_12(2)))];
end

ng = D1_ng_12 + Drug_effect;
rate = exp(ng.*12);

DATA_rate = T12_2T'*rate;

%%
% 
% save('T12_2T_v.mat','T12_2T')
% save('Forward_D24_v.mat','DATA_rate')






