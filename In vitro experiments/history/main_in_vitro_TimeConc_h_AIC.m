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
Conc_D = [0,0.125,0.25,0.5,1,2,4,8,16];
Time_D = [0,12,36];
NT_D   = length(Time_D);
NC_D   = length(Conc_D);
NR_D   = 2;
cmd  = 'CSC_DIS';



%% Load data

Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1F');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1E');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 1E');

DATA_SC = [cell2mat(Cell_SC(2:10,2:3))];
DATA_TC = [cell2mat(Cell_TC(2:10,2:3))];
DMSO_1  = cell2mat(DMSO_TC(2,6:10));
DMSO_2  = cell2mat(DMSO_TC(2,11:15));


% 
DATA_D    = zeros(NT_D,NC_D,NR_D);
% DATA(1,:,1) = mean(DMSO_1)/2.4655; % [0-48] rate = 0.0188
% DATA(1,:,2) = mean(DMSO_2)/2.4655;
% DATA(1,:,1) = mean(DMSO_1)/1.9676; % [12-48]
% DATA(1,:,2) = mean(DMSO_2)/1.9676;
% DATA(1,:,1) = mean(DMSO_1)/1.5702; % [24-48]
% DATA(1,:,2) = mean(DMSO_2)/1.5702;

% load('rate24.mat')
% DATA(1,:,1) = mean(DMSO_1)./rate;
% DATA(1,:,2) = mean(DMSO_2)./rate;

DATA_D(3,:,1) = DATA_TC(:,1)*mean(DMSO_1)/100;
DATA_D(3,:,2) = DATA_TC(:,2)*mean(DMSO_2)/100;


load('Forward_D24.mat')
DATA_D(2,:,1) = DATA_rate;
DATA_D(2,:,2) = DATA_rate;

load('T12_2T.mat')
DATA_D(1,:,1) = T12_2T;
DATA_D(1,:,2) = T12_2T;



%%  Setting the Concentration and Time points

init   = 1e4;
s      = 2;
Conc_T = [0,4,8];
Time_T = [12,24,48];
NT_T   = length(Time_T);
NC_T   = length(Conc_T);
NR_T   = 3;
cmd  = 'CSC_DIS';





%% 3 replicates

Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figures 2A');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 2B');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 2B');


DATA_SC = [cell2mat(Cell_SC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_SC(3:5,11:16))];
DATA_TC = [cell2mat(Cell_TC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_TC(3:5,11:16))];
DATA_DMSO = [cell2mat(DMSO_TC(9:11,2:9)),[nan;nan;nan],cell2mat(DMSO_TC(9:11,11:16))];
DATA_DMSO_avg = mean(DATA_DMSO,1);

DATA_T    = zeros(NT_T,NC_T,NR_T);
DATA_sc_T = zeros(NT_T,NC_T,NR_T);
for i = 1:NT_T
    k = i+2;
    DATA_T(i,:,:) = DATA_DMSO_avg(3*k-2:3*k).*DATA_TC(:,3*k-2:3*k)/100;
    DATA_sc_T(i,:,:) = squeeze(DATA_T(i,:,:)).*DATA_SC(:,3*k-2:3*k)/100;
end
% DATA(1,:,:) = init;
% DATA_sc(1,:,:) = 0.25/100*init;



r = squeeze(DATA_sc_T(1,:,:))./squeeze(DATA_T(1,:,:));

ratio = mean(r,2,'omitnan');


%% Optimization bound (DIP)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]



alpha_lb = 0;
alpha_ub = 1-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 1;
nu_lb = 0;
nu_ub = 1-1e-6;
nu_d_ub = 0;
b_beta_lb  = 1e-6;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 16;
m_lb  = 0.1;
m_ub  = 5;
c_lb  = 0;
c_ub  = 5000;

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

num_optim= 100;


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
func = @(x) get_like_alt_h(DATA_D,x,Time_D,Conc_D,NR_D,NC_D,NT_D,s,cmd) + get_like_VIP_h(DATA_T,x,Time_T,Conc_T,NR_T,NC_T,NT_T,s,cmd,ratio);
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
t = toc
[of_DIP,oi] = min(fval_hist_pe_DIP);
opt_xx_pe_DIP = params_hist_pe_DIP(oi,:);


%% Optimization bound (no DIP)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]






    
alpha_lb = 0;
alpha_ub = 1-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 1;
nu_lb = 0;
nu_ub = 1-1e-6;
nu_d_ub = 0;
b_beta_lb  = 1e-6;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 16;
m_lb  = 0.1;
m_ub  = 5;
c_lb  = 0;
c_ub  = 5000;

lb = [0,alpha_lb,beta_lb,nu_lb,1,1,1,1,1,1,...
        0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,1,1,1,c_lb];



ub = [1,alpha_ub,beta_ub,nu_ub,1,1,1,1,1,1,...
        1,alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,m_ub,1,1,1,c_ub];


A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
b = [0.1;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = [0,1];

num_optim= 100;


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
func = @(x) get_like_alt_h(DATA_D,x,Time_D,Conc_D,NR_D,NC_D,NT_D,s,cmd) + get_like_VIP_h(DATA_T,x,Time_T,Conc_T,NR_T,NC_T,NT_T,s,cmd,ratio) ;
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
t = toc
[of_nDIP,oi] = min(fval_hist_pe_nDIP);
opt_xx_pe_nDIP = params_hist_pe_nDIP(oi,:);




%% Compute the AIC


nDIP_AIC = 2*9 + 2*of_nDIP;
DIP_AIC = 2*12 + 2*of_DIP;



    

save_name = strcat('Result/In_vitro_TimeConc_122448_h_AIC.mat');

save(save_name)


