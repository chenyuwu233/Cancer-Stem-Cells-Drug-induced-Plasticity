%%  AIC script for COLO858--Vemurafenib:
% - DIP: Drug-induced plasticity
% - NP: Natural
% - SC: Stem-like cells Cytotoxic
% 
parpool('local',50)
warning('off','MATLAB:integral:NonFiniteValue')
% addpath('C:\Users\euclid\OneDrive\UMN-My-gear\Phenotypic_Switching_model\Switching_model(MATLAB)\Main functions')

%%
% seed_num = 52;
% 
% rng(seed_num)



load('DATA_COLO858_60h.mat')

%%  Setting the Concentration and Time points


DATA_full = DATA_COLO858_ss;

Conc_full = Conc;
% ss_idx = [1,2,4,6]; % Include the DMSO
% Conc = Conc(ss_idx);
DATA = DATA_full;


s      = 2;
[NT,NC,NR] = size(DATA);
cmd  = 'CSC_DIS';




%%








%% Optimization bound (Hill2_switching_death)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]


alpha_lb = 0;
alpha_ub = 0.1-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 0.1;
nu_lb = 0;
nu_ub = 0.1;
nu_d_ub = 0.1;
b_beta_lb  = 1e-6;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 50;
tv1_lb = 0;
tv1_ub = 5;
tv_lb = 0;
tv_ub = 36;
c_lb  = 0;
c_ub  = 50;

lb_DIP_Asy = [0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,1,1,...
        1,alpha_lb,beta_lb,0,b_beta_lb,E_lb,b_nu_lb,E_lb,tv1_lb,tv_lb,c_lb];



ub_DIP_Asy = [0,alpha_ub,beta_ub,0.5,b_beta_ub,E_ub,1,1,...
        1,alpha_ub,beta_ub,0.5,b_beta_ub,E_ub,b_nu_ub,E_ub,tv1_ub,tv_ub,c_ub];

% 
% A = [0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0];
% b = [0;0];

Aeq = [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = 1;

num_optim= 50;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub_DIP_Asy)).*(ub_DIP_Asy-lb_DIP_Asy)+lb_DIP_Asy;
    x_init = [x_init;xi];
end







%% Optimization (Point estimate DIP)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_TD_meanflow(DATA,x,Time_ss,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_DIP_Asy   = [];
params_hist_pe_DIP_Asy = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),[],[],[],[],lb_DIP,ub_DIP,[],options1); 
    fval_hist_pe_DIP_Asy = [fval_hist_pe_DIP_Asy,ff];
    params_hist_pe_DIP_Asy = [params_hist_pe_DIP_Asy;xx];
%     catch
%     end
end
[of_DIP_Asy,oi] = min(fval_hist_pe_DIP_Asy);
opt_xx_pe_DIP_Asy = params_hist_pe_DIP_Asy(oi,:);


%%


lb_DIP_nAsy = [0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,1,1,...
        1,alpha_lb,beta_lb,0,b_beta_lb,E_lb,b_nu_lb,E_lb,tv1_lb,tv_lb,c_lb];



ub_DIP_nAsy = [0,alpha_ub,beta_ub,0,b_beta_ub,E_ub,1,1,...
        1,alpha_ub,beta_ub,0,b_beta_ub,E_ub,b_nu_ub,E_ub,tv1_ub,tv_ub,c_ub];

x_init_DIP_nAsy = x_init;
x_init_DIP_nAsy(:,4) = 0;
x_init_DIP_nAsy(:,12) = 0;

%% Optimization (Point estimate nDIP)


options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
% func = @(x) get_like_TD_alt2(DATA,x,Time_ss,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_DIP_nAsy   = [];
params_hist_pe_DIP_nAsy = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_DIP_nAsy(j,:),[],[],[],[],lb_DIP_nNP,ub_DIP_nNP,[],options1); 
    fval_hist_pe_DIP_nAsy = [fval_hist_pe_DIP_nAsy,ff];
    params_hist_pe_DIP_nAsy = [params_hist_pe_DIP_nAsy;xx];
%     catch
%     end
end
[of_DIP_nAsy,oi] = min(fval_hist_pe_DIP_nAsy);
opt_xx_pe_DIP_nAsy = params_hist_pe_DIP_nAsy(oi,:);



%%


lb_nDIP_Asy = [0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,1,1,...
        1,alpha_lb,beta_lb,0,b_beta_lb,E_lb,1,1,tv1_lb,tv_lb,c_lb];



ub_nDIP_Asy = [0,alpha_ub,beta_ub,0.5,b_beta_ub,E_ub,1,1,...
        1,alpha_ub,beta_ub,0.5,b_beta_ub,E_ub,1,1,tv1_ub,tv_ub,c_ub];

x_init_nDIP_Asy = x_init;
x_init_nDIP_Asy(:,15) = 1;
x_init_nDIP_Asy(:,16) = 1;

%% Optimization (Point estimate nDIP)


options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
% func = @(x) get_like_TD_alt2(DATA,x,Time_ss,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_nDIP_Asy   = [];
params_hist_pe_nDIP_Asy = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_nDIP_Asy(j,:),[],[],[],[],lb_nDIP,ub_nDIP,[],options1); 
    fval_hist_pe_nDIP_Asy = [fval_hist_pe_nDIP_Asy,ff];
    params_hist_pe_nDIP_Asy = [params_hist_pe_nDIP_Asy;xx];
%     catch
%     end
end
[of_nDIP_Asy,oi] = min(fval_hist_pe_nDIP_Asy);
opt_xx_pe_nDIP_Asy = params_hist_pe_nDIP_Asy(oi,:);


%%


lb_nDIP_nAsy = [0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,1,1,...
        1,alpha_lb,beta_lb,0,b_beta_lb,E_lb,1,1,tv1_lb,tv_lb,c_lb];



ub_nDIP_nAsy = [0,alpha_ub,beta_ub,0,b_beta_ub,E_ub,1,1,...
        1,alpha_ub,beta_ub,0,b_beta_ub,E_ub,1,1,tv1_ub,tv_ub,c_ub];

x_init_nDIP_nAsy = x_init;
x_init_nDIP_nAsy(:,15) = 1;
x_init_nDIP_nAsy(:,16) = 1;
x_init_nDIP_nAsy(:,4)  = 0;
x_init_nDIP_nAsy(:,12) = 0;

%% Optimization (Point estimate nDIP)


options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
% func = @(x) get_like_TD_alt2(DATA,x,Time_ss,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_nDIP_nAsy   = [];
params_hist_pe_nDIP_nAsy = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_nDIP_nAsy(j,:),[],[],[],[],lb_nDIP_nAsy,ub_nDIP_nAsy,[],options1); 
    fval_hist_pe_nDIP_nAsy = [fval_hist_pe_nDIP_nAsy,ff];
    params_hist_pe_nDIP_nAsy = [params_hist_pe_nDIP_nAsy;xx];
%     catch
%     end
end
[of_nDIP_nNP,oi] = min(fval_hist_pe_nDIP_nAsy);
opt_xx_pe_nDIP_nNP = params_hist_pe_nDIP_nAsy(oi,:);





%% Compute the AIC:

num_par_DIP_Asy = sum(lb_DIP_Asy ~= ub_DIP_Asy);
num_par_DIP_nAsy = sum(lb_DIP_nAsy ~= ub_DIP_nAsy);
num_par_nDIP_Asy = sum(lb_nDIP_Asy ~= ub_nDIP_Asy);
num_par_nDIP_nAsy = 8;

AIC_DIP = 2*num_par_DIP_Asy + 2*of_DIP;
AIC_DIP_nNP = 2*num_par_DIP_nAsy + 2*of_DIP_nNP;
AIC_nDIP = 2*num_par_nDIP_Asy + 2*of_nDIP; 
AIC_nDIP_nNP = 2*num_par_nDIP_nAsy + 2*of_nDIP_nNP;

AICc_DIP = AIC_DIP + (2*num_par_DIP_Asy^2 + 2*num_par_DIP_Asy)/(29*6*4 - num_par_DIP_Asy - 1);
AICc_DIP_nNP = AIC_DIP_nNP + (2*num_par_DIP_nAsy^2 + 2*num_par_DIP_nAsy)/(29*6*4 - num_par_DIP_nAsy - 1);
AICc_nDIP = AIC_nDIP + (2*num_par_nDIP_Asy^2 + 2*num_par_nDIP_Asy)/(29*6*4 - num_par_nDIP_Asy - 1);
AICc_nDIP_nNP = AIC_nDIP_nNP + (2*num_par_nDIP_nAsy^2 + 2*num_par_nDIP_nAsy)/(29*6*4 - num_par_nDIP_nAsy - 1);


    

save_name = strcat('Result/In_vitro_COLO858_60h_AIC.mat');

save(save_name)


