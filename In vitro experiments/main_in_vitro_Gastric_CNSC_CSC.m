%%  Generating the theta:
% - s: sub-population
% - p: Initial proportion
% - 
% % 


parpool('local',20)
warning('off','MATLAB:integral:NonFiniteValue')



%%  Setting the Concentration and Time points

% init   = 5000;
s      = 2;
Conc = [0,0.125,0.25,0.5,1,2,4,8,16];
Time = [0,12,36];
NT   = length(Time);
NC   = length(Conc);
NR   = 2;
cmd  = 'CSC_DIS';
lam_s_range       = [0,0.1];
lam_d_range       = [0,0.5];



%% Load data

Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1F');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1E');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 1E');

DATA_SC = [cell2mat(Cell_SC(2:10,2:3))];
DATA_TC = [cell2mat(Cell_TC(2:10,2:3))];
DMSO_1  = cell2mat(DMSO_TC(2,6:10));
DMSO_2  = cell2mat(DMSO_TC(2,11:15));


% 
DATA    = zeros(NT,NC,NR);
% DATA(1,:,1) = mean(DMSO_1)/2.4655; % [0-48] rate = 0.0188
% DATA(1,:,2) = mean(DMSO_2)/2.4655;
% DATA(1,:,1) = mean(DMSO_1)/1.9676; % [12-48]
% DATA(1,:,2) = mean(DMSO_2)/1.9676;
% DATA(1,:,1) = mean(DMSO_1)/1.5702; % [24-48]
% DATA(1,:,2) = mean(DMSO_2)/1.5702;

% load('rate24.mat')
% DATA(1,:,1) = mean(DMSO_1)./rate;
% DATA(1,:,2) = mean(DMSO_2)./rate;

DATA(3,:,1) = DATA_TC(:,1)*mean(DMSO_1)/100;
DATA(3,:,2) = DATA_TC(:,2)*mean(DMSO_2)/100;


load('Forward_D24_v.mat')
DATA(2,:,1) = DATA_rate(1,:);
DATA(2,:,2) = DATA_rate(2,:);

load('T12_2T_v.mat')
DATA(1,:,1) = T12_2T(1);
DATA(1,:,2) = T12_2T(2);

DATA = round(DATA);



%% Optimization bound (DIP)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]



alpha_lb = 0;
alpha_ub = 1-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 1;
nu_lb = 0;
nu_ub = 1-1e-6;
nu_d_ub = 1-1e-6;
b_beta_lb  = 1e-6;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 16;
m_lb  = 0.1;
m_ub  = 5;
c_lb  = 0;
c_ub  = 300;

lb = [0,alpha_lb,beta_lb,nu_lb,1,1,1,1,1,1,...
        0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,b_nu_lb,E_lb,m_lb,c_lb];



ub = [1,alpha_ub,beta_ub,nu_ub,1,1,1,1,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,m_ub,b_nu_ub,E_ub,m_ub,c_ub];


A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
b = [0.1;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = [0,1];

num_optim= 20;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(2) = xi(3) + rand*0.1;
    xi(12) = xi(13);
    xi(1)  = 0;
    xi(11) = 1;
    x_init = [x_init;xi];
end
    




%% Optimization (Point estimate death)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
func = @(x) get_like_alt_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_DIP   = [];
params_hist_pe_DIP = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_DIP = [fval_hist_pe_DIP,ff];
    params_hist_pe_DIP = [params_hist_pe_DIP;xx];
%     catch
%     end
end
t_DIP = toc
[of_DIP,oi] = min(fval_hist_pe_DIP);
opt_xx_pe_DIP = params_hist_pe_DIP(oi,:);

%% Optimization bound (DIP) TD2
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]



alpha_lb = 0;
alpha_ub = 1-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 1;
nu_lb = 0;
nu_ub = 1-1e-6;
nu_d_ub = 1-1e-6;
b_beta_lb  = 1e-6;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 16;
m_lb  = 0.1;
m_ub  = 5;
tv1_lb = 0;
tv1_ub = 1;
tv_lb = 0;
tv_ub = 36;
c_lb  = 0;
c_ub  = 300;

lb = [0,alpha_lb,beta_lb,nu_lb,1,1,1,1,1,1,...
        0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,b_nu_lb,E_lb,m_lb,tv1_lb,tv_lb,c_lb];



ub = [1,alpha_ub,beta_ub,nu_ub,1,1,1,1,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,m_ub,b_nu_ub,E_ub,m_ub,tv1_ub,tv_ub,c_ub];


A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
b = [0.1;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0];
beq = [0,1];

num_optim= 20;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(2) = xi(3) + rand*0.1;
    xi(12) = xi(13);
    xi(1)  = 0;
    xi(11) = 1;
    x_init = [x_init;xi];
end
    




%% Optimization (Point estimate death)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
func = @(x) get_like_TD_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_TD2   = [];
params_hist_pe_TD2 = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_TD2 = [fval_hist_pe_TD2,ff];
    params_hist_pe_TD2 = [params_hist_pe_TD2;xx];
%     catch
%     end
end
t_TD2 = toc
[of_TD2,oi] = min(fval_hist_pe_TD2);
opt_xx_pe_TD2 = params_hist_pe_TD2(oi,:);

%% Optimization bound (no DIP)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]






    
alpha_lb = 0;
alpha_ub = 1-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 1;
nu_lb = 0;
nu_ub = 1-1e-6;
nu_d_ub = 1-1e-6;
b_beta_lb  = 1e-6;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 16;
m_lb  = 0.1;
m_ub  = 5;
c_lb  = 0;
c_ub  = 300;

lb = [0,1,1,1,1,1,1,1,1,1,...
        0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,1,1,1,c_lb];



ub = [1,1,1,1,1,1,1,1,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,m_ub,1,1,1,c_ub];


A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
b = [0.1;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = [0,1];

num_optim= 20;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(2) = xi(3) + rand*0.1;
    xi(12) = xi(13);
    xi(1)  = 0;
    xi(11) = 1;
    x_init = [x_init;xi];
end
    






%% Optimization (Point estimate death)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
func = @(x) get_like_alt_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_nDIP   = [];
params_hist_pe_nDIP = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_nDIP = [fval_hist_pe_nDIP,ff];
    params_hist_pe_nDIP = [params_hist_pe_nDIP;xx];
%     catch
%     end
end
t_nDIP = toc
[of_nDIP,oi] = min(fval_hist_pe_nDIP);
opt_xx_pe_nDIP = params_hist_pe_nDIP(oi,:);


%% Logistic time delay with 1 parameter







    
alpha_lb = 0;
alpha_ub = 1-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 1;
nu_lb = 0;
nu_ub = 1-1e-6;
nu_d_ub = 1-1e-6;
b_beta_lb  = 1e-6;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 16;
m_lb  = 0.1;
m_ub  = 5;
tv1_lb = 0;
tv1_ub = 1;
tv_lb = 0;
tv_ub = 36;
c_lb  = 0;
c_ub  = 300;

lb = [0,1,1,1,1,1,1,1,1,1,...
        0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,1,1,1,tv1_lb,c_lb];



ub = [1,1,1,1,1,1,1,1,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,m_ub,1,1,1,tv1_ub,c_ub];


A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
b = [0.1;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
beq = [0,1];

num_optim= 20;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(2) = xi(3) + rand*0.1;
    xi(12) = xi(13);
    xi(1)  = 0;
    xi(11) = 1;
    x_init = [x_init;xi];
end
    






%% Optimization (Point estimate death)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
func = @(x) get_like_TD_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_nTD1   = [];
params_hist_pe_nTD1 = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_nTD1 = [fval_hist_pe_nTD1,ff];
    params_hist_pe_nTD1 = [params_hist_pe_nTD1;xx];
%     catch
%     end
end
t_nTD1 = toc
[of_nTD1,oi] = min(fval_hist_pe_nTD1);
opt_xx_pe_nTD1 = params_hist_pe_nTD1(oi,:);


%% Logistic time delay with 2 parameters 







    
alpha_lb = 0;
alpha_ub = 1-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 1;
nu_lb = 0;
nu_ub = 1-1e-6;
nu_d_ub = 1-1e-6;
b_beta_lb  = 1e-6;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 16;
m_lb  = 0.1;
m_ub  = 5;
tv1_lb = 0;
tv1_ub = 1;
tv_lb = 0;
tv_ub = 36;
c_lb  = 0;
c_ub  = 300;

lb = [0,1,1,1,1,1,1,1,1,1,...
        0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,1,1,1,tv1_lb,tv_lb,c_lb];



ub = [1,1,1,1,1,1,1,1,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,m_ub,1,1,1,tv1_ub,tv_ub,c_ub];


A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
b = [0.1;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0];
beq = [0,1];

num_optim= 20;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(2) = xi(3) + rand*0.1;
    xi(12) = xi(13);
    xi(1)  = 0;
    xi(11) = 1;
    x_init = [x_init;xi];
end
    






%% Optimization (Point estimate death)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
func = @(x) get_like_TD_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_nTD2   = [];
params_hist_pe_nTD2 = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_nTD2 = [fval_hist_pe_nTD2,ff];
    params_hist_pe_nTD2 = [params_hist_pe_nTD2;xx];
%     catch
%     end
end
t_nTD2 = toc
[of_nTD2,oi] = min(fval_hist_pe_nTD2);
opt_xx_pe_nTD2 = params_hist_pe_nTD2(oi,:);


%% Compute the AIC

DIP_num_par = 12;
nDIP_num_par = 6;
nDIP_TD1_num_par = 7;
nDIP_TD2_num_par = 8;
DIP_TD2_num_par = 14;

nDIP_AIC = 2*6 + 2*of_nDIP;
DIP_AIC = 2*DIP_num_par + 2*of_DIP;
nDIP_TD1_AIC = 2*7 + 2*of_nTD1;
nDIP_TD2_AIC = 2*8 + 2*of_nTD2;
DIP_TD2_AIC = 2*14 + 2*of_TD2;


DIP_AICc = DIP_AIC + (2*DIP_num_par^2 + 2*DIP_num_par)/(3*9*2 - DIP_num_par - 1);
nDIP_AICc = nDIP_AIC + (2*nDIP_num_par^2 + 2*nDIP_num_par)/(3*9*2 - nDIP_num_par - 1);
nDIP_TD1_AICc = nDIP_TD1_AIC + (2*nDIP_TD1_num_par^2 + 2*nDIP_TD1_num_par)/(3*9*2 - nDIP_TD1_num_par - 1);
nDIP_TD2_AICc = nDIP_TD2_AIC + (2*nDIP_TD2_num_par^2 + 2*nDIP_TD2_num_par)/(3*9*2 - nDIP_TD2_num_par - 1);
DIP_TD2_AICc = DIP_TD2_AIC + (2*DIP_TD2_num_par^2 + 2*DIP_TD2_num_par)/(3*9*2 - DIP_TD2_num_par - 1);




    

save_name = strcat('Result/In_vitro_Conc_122448_TD_h_AIC.mat');

save(save_name)


