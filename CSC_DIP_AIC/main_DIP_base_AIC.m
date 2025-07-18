%%  Generating the theta:
% - s: sub-population
% - p: Initial proportion
% - 
% 
parpool('local',20)
warning('off','MATLAB:integral:NonFiniteValue')
% addpath('C:\Users\euclid\OneDrive\UMN-My-gear\Phenotypic_Switching_model\Switching_model(MATLAB)\Main functions')

%%
seed_num = 26;

rng(seed_num)





%%  Setting the Concentration and Time points

init   = 1000;
s      = 2;
Conc = 10^(6)*[0  31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6)];
NT   = 13;
Time = [0 : NT-1]*3;
NC   = length(Conc);
NR   = 20;
cmd  = 'CSC_DIS';

% p_range      = [0.25,0.75];
beta_s_range      = [1e-3,0.9];
beta_d_range      = [1e-3,0.5];
lam_s_range       = [0,0.1];
lam_d_range       = [0,0.05];
nu_dediff_range = [0,0.03];
b_beta_range    = [0.8,0.9];
b_nu_range      = [1+1e-6,1.1];
E_range         = [0.0625,2.5];
c_range         = [0,10];


c      = rand*10;
beta1  = rand*(beta_s_range(2)-beta_s_range(1)) + beta_s_range(1);
beta2  = rand*(beta_d_range(2)-beta_d_range(1)) + beta_d_range(1);
alpha1 = beta1 + rand*lam_s_range(2);
alpha2 = beta2 + rand*lam_d_range(2);
nu12   = rand*nu_dediff_range(2);
nu21   = 0;
b1_temp = rand*(b_beta_range(2)-b_beta_range(1)) + b_beta_range(1);
b2_temp = rand*(b_beta_range(2)-b_beta_range(1)) + b_beta_range(1);
b1_beta = max(b1_temp,b2_temp);
b2_beta = min(b1_temp,b2_temp);
% b2_nu   = rand*(b_nu_range(2)-b_nu_range(1)) + b_nu_range(1);
b2_nu   = 1;
E1_beta = rand*(E_range(2)-E_range(1)) + E_range(1);
E2_beta = rand*(E_range(2)-E_range(1)) + E_range(1);
E2_nu = rand*(E_range(2)-E_range(1)) + E_range(1);

p_1     = 0;
theta   = [p_1,alpha1,beta1,nu12,b1_beta,E1_beta,1,1,...
              1-p_1,alpha2,beta2,nu21,b2_beta,E2_beta,b2_nu,E2_nu,c]

A = [alpha1-beta1,nu12;nu21,alpha2-beta2];


eig(A)
sp = get_stable_p(A);
p_diff = p_1 -sp(1);




%%


tic
DATA = Switching_gen_ip(init,theta,Time,Conc,NR,NC,NT,s,cmd);
t_gen = toc

tic
opt_fval = get_like_ip(DATA,theta,Time,Conc,NR,NC,NT,s,cmd);
t_like = toc









alpha_lb = 0;
alpha_ub = 1-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 1;
nu_lb = 0;
nu_ub = 1-1e-6;
nu_d_ub = 0;
b_beta_lb  = 0.5;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 1.5;
E_lb  = 1e-6;
E_ub  = 5;
c_lb  = 0;
c_ub  = 10;


%% Optimization (DIP_Asy)
lb_DIP_Asy = [0,alpha_lb,beta_lb,nu_lb,b_beta_lb,E_lb,1,1,...
        1,alpha_lb,beta_lb,nu_lb,b_beta_lb,E_lb,b_nu_lb,E_lb,c_lb];



ub_DIP_Asy = [0,alpha_ub,beta_ub,nu_ub,b_beta_ub,E_ub,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,b_nu_ub,E_ub,c_ub];


% A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%      0,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,0,0];
% b = [lam_s_range(2);lam_d_range(2)];

Aeq = [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0];
beq = 1;

num_optim= 20;


x_init_DIP_Asy = [];
for i = 1:num_optim
    xi = rand(1,length(ub_DIP_Asy)).*(ub_DIP_Asy-lb_DIP_Asy)+lb_DIP_Asy;
    xi(2) = xi(3) + rand*lam_s_range(2);
    xi(4) = rand*lam_d_range(2);
    xi(10) = xi(11);
    x_init_DIP_Asy = [x_init_DIP_Asy;xi];
end




options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_ip(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_DIP_Asy(j,:),[],[],Aeq,beq,lb_DIP_Asy,ub_DIP_Asy,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
%     catch
%     end
end
[of_DIP_Asy,oi] = min(fval_hist_pe);
opt_xx_pe_DIP_Asy = params_hist_pe(oi,:);



ub_DIP_Asy(end) = 1e4;

func = @(x) get_like_LLM_ip(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_DIP_Asy(j,:),[],[],Aeq,beq,lb_DIP_Asy,ub_DIP_Asy,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
%     catch
%     end
end
[of_LLM_DIP_Asy,oi] = min(fval_hist_pe);
opt_xx_pe_LLM_DIP_Asy = params_hist_pe(oi,:);

%% Optimization (DIP_nAsy)

lb_DIP_nAsy = [0,alpha_lb,beta_lb,nu_lb,b_beta_lb,E_lb,1,1,...
        1,alpha_lb,beta_lb,nu_lb,b_beta_lb,E_lb,b_nu_lb,E_lb,c_lb];



ub_DIP_nAsy = [0,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,b_nu_ub,E_ub,c_ub];


% A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%      0,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,0,0];
% b = [lam_s_range(2);lam_d_range(2)];

Aeq = [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0];
beq = 1;

num_optim= 20;


x_init_DIP_nAsy = [];
for i = 1:num_optim
    xi = rand(1,length(ub_DIP_nAsy)).*(ub_DIP_nAsy-lb_DIP_nAsy)+lb_DIP_nAsy;
    xi(2) = xi(3) + rand*lam_s_range(2);
    xi(4) = rand*lam_d_range(2);
    xi(10) = xi(11);
    x_init_DIP_nAsy = [x_init_DIP_nAsy;xi];
end




options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_ip(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_DIP_nAsy(j,:),[],[],Aeq,beq,lb_DIP_nAsy,ub_DIP_nAsy,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
%     catch
%     end
end
[of_DIP_nAsy,oi] = min(fval_hist_pe);
opt_xx_pe_DIP_nAsy = params_hist_pe(oi,:);

ub_DIP_nAsy(end) = 1e4;

func = @(x) get_like_LLM_ip(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_DIP_nAsy(j,:),[],[],Aeq,beq,lb_DIP_nAsy,ub_DIP_nAsy,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
%     catch
%     end
end
[of_LLM_DIP_nAsy,oi] = min(fval_hist_pe);
opt_xx_pe_LLM_DIP_nAsy = params_hist_pe(oi,:);


%% Optimization (nDIP_Asy)

lb_nDIP_Asy = [0,alpha_lb,beta_lb,nu_lb,b_beta_lb,E_lb,1,1,...
        1,alpha_lb,beta_lb,nu_lb,b_beta_lb,E_lb,1,1,c_lb];



ub_nDIP_Asy = [0,alpha_ub,beta_ub,nu_ub,b_beta_ub,E_ub,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,1,1,c_ub];


% A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%      0,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,0,0];
% b = [lam_s_range(2);lam_d_range(2)];

Aeq = [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0];
beq = 1;

num_optim= 20;


x_init_nDIP_Asy = [];
for i = 1:num_optim
    xi = rand(1,length(ub_nDIP_Asy)).*(ub_nDIP_Asy-lb_nDIP_Asy)+lb_nDIP_Asy;
    xi(2) = xi(3) + rand*lam_s_range(2);
    xi(4) = rand*lam_d_range(2);
    xi(10) = xi(11);
    x_init_nDIP_Asy = [x_init_nDIP_Asy;xi];
end




options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_ip(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_nDIP_Asy(j,:),[],[],Aeq,beq,lb_nDIP_Asy,ub_nDIP_Asy,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
%     catch
%     end
end
[of_nDIP_Asy,oi] = min(fval_hist_pe);
opt_xx_pe_nDIP_Asy = params_hist_pe(oi,:);


ub_nDIP_Asy(end) = 1e4;

func = @(x) get_like_LLM_ip(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_nDIP_Asy(j,:),[],[],Aeq,beq,lb_nDIP_Asy,ub_nDIP_Asy,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
%     catch
%     end
end
[of_LLM_nDIP_Asy,oi] = min(fval_hist_pe);
opt_xx_pe_LLM_nDIP_Asy = params_hist_pe(oi,:);
   


%% Optimization (nDIP_nAsy)

lb_nDIP_nAsy = [0,alpha_lb,beta_lb,nu_lb,b_beta_lb,E_lb,1,1,...
        1,alpha_lb,beta_lb,nu_lb,b_beta_lb,E_lb,1,1,c_lb];



ub_nDIP_nAsy = [0,alpha_lb,beta_lb,nu_d_ub,b_beta_lb,E_lb,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,1,1,c_ub];


% A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%      0,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,0,0];
% b = [lam_s_range(2);lam_d_range(2)];

Aeq = [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0];
beq = 1;

num_optim= 20;


x_init_nDIP_nAsy = [];
for i = 1:num_optim
    xi = rand(1,length(ub_nDIP_nAsy)).*(ub_nDIP_nAsy-lb_nDIP_nAsy)+lb_nDIP_nAsy;
    xi(2) = xi(3) + rand*lam_s_range(2);
    xi(4) = rand*lam_d_range(2);
    xi(10) = xi(11);
    x_init_nDIP_nAsy = [x_init_nDIP_nAsy;xi];
end





options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_ip(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_nDIP_nAsy(j,:),[],[],Aeq,beq,lb_nDIP_nAsy,ub_nDIP_nAsy,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
%     catch
%     end
end
[of_nDIP_nAsy,oi] = min(fval_hist_pe);
opt_xx_pe_nDIP_nAsy = params_hist_pe(oi,:);


ub_nDIP_nAsy(end) = 1e4;

func = @(x) get_like_LLM_ip(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];

parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init_nDIP_nAsy(j,:),[],[],Aeq,beq,lb_nDIP_nAsy,ub_nDIP_nAsy,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
%     catch
%     end
end
[of_LLM_nDIP_nAsy,oi] = min(fval_hist_pe);
opt_xx_pe_LLM_nDIP_nAsy = params_hist_pe(oi,:);


%%

of_DIP_Asy = func(opt_xx_pe_DIP_Asy);


num_par_DIP_Asy = sum(lb_DIP_Asy ~= ub_DIP_Asy);
num_par_DIP_nAsy = sum(lb_DIP_nAsy ~= ub_DIP_nAsy);
num_par_nDIP_Asy = 5;
num_par_nDIP_nAsy = 5;
num_par_LLM_DIP_Asy = num_par_DIP_Asy - 2;
num_par_LLM_DIP_nAsy = num_par_DIP_nAsy - 2;
num_par_LLM_nDIP_Asy = num_par_nDIP_Asy - 1;
num_par_LLM_nDIP_nAsy = num_par_nDIP_nAsy - 1;


AIC_DIP_Asy = 2*num_par_DIP_Asy + 2*of_DIP_Asy;
AIC_DIP_nAsy = 2*num_par_DIP_nAsy + 2*of_DIP_nAsy;
AIC_nDIP_Asy = 2*num_par_nDIP_Asy + 2*of_nDIP_Asy; 
AIC_nDIP_nAsy = 2*num_par_nDIP_nAsy + 2*of_nDIP_nAsy;
AIC_LLM_DIP_Asy = 2*num_par_LLM_DIP_Asy + 2*of_LLM_DIP_Asy;
AIC_LLM_DIP_nAsy = 2*num_par_LLM_DIP_nAsy + 2*of_LLM_DIP_nAsy;
AIC_LLM_nDIP_Asy = 2*num_par_LLM_nDIP_Asy + 2*of_LLM_nDIP_Asy; 
AIC_LLM_nDIP_nAsy = 2*num_par_LLM_nDIP_nAsy + 2*of_LLM_nDIP_nAsy;


AICc_DIP_Asy = AIC_DIP_Asy + (2*num_par_DIP_Asy^2 + 2*num_par_DIP_Asy)/((NT-1)*NC*NR - num_par_DIP_Asy - 1);
AICc_DIP_nAsy = AIC_DIP_nAsy + (2*num_par_DIP_nAsy^2 + 2*num_par_DIP_nAsy)/((NT-1)*NC*NR - num_par_DIP_nAsy - 1);
AICc_nDIP_Asy = AIC_nDIP_Asy + (2*num_par_nDIP_Asy^2 + 2*num_par_nDIP_Asy)/((NT-1)*NC*NR - num_par_nDIP_Asy - 1);
AICc_nDIP_nAsy = AIC_nDIP_nAsy + (2*num_par_nDIP_nAsy^2 + 2*num_par_nDIP_nAsy)/((NT-1)*NC*NR - num_par_nDIP_nAsy - 1);


save_name = strcat('Result/nDIP_base_LLM_AIC_',num2str(seed_num),'.mat');

save(save_name)


