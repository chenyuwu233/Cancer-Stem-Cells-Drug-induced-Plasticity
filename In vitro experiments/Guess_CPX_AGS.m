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

DATA    = zeros(NT,NC,NR);
DATA_sc = zeros(NT,NC,NR);
for i = 1:5
    k = i;
    DATA(i,:,:) = DATA_DMSO_avg(3*k-2:3*k).*DATA_TC(:,3*k-2:3*k)/100;
    DATA_sc(i,:,:) = squeeze(DATA(i,:,:)).*DATA_SC(:,3*k-2:3*k)/100;
end

E_data = mean(DATA(:,1,:),3,'omitnan');

Fit_5 = fit(Time',E_data,'exp1');

%% data From Alvaro resource

Conc = [0.002,0.0041,0.0081,0.016,0.032,0.065,0.13,0.26,0.52,2.1,8.3,17,33,66];
Inhi = [5.55,6.37,2.09,3,6.78,2.93,3.71,5.54,5.07,1.18,51.38,66.25,76.39,92.94];
% Conc = [0.0081,8.3,17,33,66];
% Inhi = [2.09,51.38,66.25,76.39,92.94];


lb = [0,0,1e-6,1];
ub = [10,10,10,1];
num_optim = 20;
x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    x_init = [x_init;xi];
end

ub = [10,100,66,1];




func = @(x) Hill4_fit(x,Conc,Inhi);

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','iter','algorithm','sqp');

fval_hist_pe = [];
params_hist_pe = [];

for j = 1:num_optim

    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),[],[],[],[],lb,ub,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];

end

[of,oi] = min(fval_hist_pe);
opt_xx_pe = params_hist_pe(oi,:)


% opt_xx_pe_hist = [];
% for i = 1:100
%     idx = randsample(14,10,1);
%     Conc_i = [];
%     Inhi_i = [];
%     for j = 1:10
%         Conc_i = [Conc_i,Conc(idx(j))];
%         Inhi_i = [Inhi_i,Inhi(idx(j))];
%     end
%     func = @(x) Hill4_fit(x,Conc_i,Inhi_i);
% 
%     fval_hist_pe = [];
%     params_hist_pe = [];
% 
%     for j = 1:num_optim
% 
%         [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),[],[],[],[],lb,ub,[],options1); 
%         fval_hist_pe = [fval_hist_pe,ff];
%         params_hist_pe = [params_hist_pe;xx];
% 
%     end
% 
%     [of,oi] = min(fval_hist_pe);
%     ope = params_hist_pe(oi,:);
%     opt_xx_pe_hist = [opt_xx_pe_hist;ope];
% 
% end
% est_hist = [est_hist;opt_xx_pe];

%%
CI_fit = confint(Fit_5);
CI_b = CI_fit(:,2);

br_vec = [CI_b(1),Fit_5.b,CI_b(2)];
est_hist = [];
for gbr = br_vec
% gbr = Fit_5.b;
% gbr = CI_b(1);
% gbr = CI_b(2);
    br  = gbr.*(100-Inhi)./100;




% br = [init, gbr*(100-51.38)/100, gbr*(100-66.25)/100, gbr*(100-76.39)/100, gbr*(100-92.937)/100];

    Drug_effect = br-gbr;

%%
% %% data From Paper
% 
% Conc = [0,0.125,0.25,0.5,1,2];
% 
% Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1F');
% Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1E');
% DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 1E');
% 
% DATA_SC = [cell2mat(Cell_SC(2:10,2:3))];
% DATA_TC = [cell2mat(Cell_TC(2:10,2:3))];
% DMSO_1  = cell2mat(DMSO_TC(2,6:10));
% DMSO_2  = cell2mat(DMSO_TC(2,11:15));
% 
% DATA    = zeros(2,9,2);
% % DATA(1,:,1) = mean(DMSO_1)/4.2931;
% % DATA(1,:,2) = mean(DMSO_2)/4.2931;
% DATA(1,:,1) = mean(DMSO_1)/1.6096; % [12-48]
% DATA(1,:,2) = mean(DMSO_2)/1.6096;
% 
% 
% DATA(2,:,1) = DATA_TC(:,1)*mean(DMSO_1)/100;
% DATA(2,:,2) = DATA_TC(:,2)*mean(DMSO_2)/100;
% 
% gr_est = [];
% for i = 1:length(Conc)
%     gr_est = [gr_est,log(mean(DATA(2,i,:))/mean(DATA(1,i,:)))/36];
% end
% 
% Drug_effect = gr_est - gr_est(1);
%% Fit 2 par Hill x = [b,E]


func = @(x) Hill_fit(x,Conc,Drug_effect);

lb = [0.5,1e-6,1e-6];
ub = [1,60,5];
num_optim = 20;
x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    x_init = [x_init;xi];
end

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','iter','algorithm','sqp');

fval_hist_pe = [];
params_hist_pe = [];

for j = 1:num_optim

    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),[],[],[],[],lb,ub,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];

end

[of,oi] = min(fval_hist_pe);
opt_xx_pe = params_hist_pe(oi,:);
est_hist = [est_hist;opt_xx_pe];

end

% %% Plot the fit
% 
% 
% X = linspace(0,70,140);
% Y = zeros(1,140);
% for i = 1:140
%     Y(i) = log(get_Hill2(X(i),opt_xx_pe(1),opt_xx_pe(2)));
% end
% plot(X,Y)
% hold on
% scatter(Conc,Drug_effect)
% xlabel('Concentration')
% ylabel('Drug effect on net growth rate')
