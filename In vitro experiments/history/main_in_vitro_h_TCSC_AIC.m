%%  Generating the theta:
% - s: sub-population
% - p: Initial proportion
% - 
% % 


parpool('local',50)
warning('off','MATLAB:integral:NonFiniteValue')



%%  Setting the Concentration and Time points

init   = 1e4;
s      = 2;
Conc = [0,4,8];
Time = [2,6,12,24,48];
NT   = length(Time);
NC   = length(Conc);
NR   = 3;
cmd  = 'CSC_DIS';
lam_s_range       = [0,0.1];
lam_d_range       = [0,0.5];



%% Load data

%% 3 replicates

Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figures 2A');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 2B');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 2B');


DATA_SC = [cell2mat(Cell_SC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_SC(3:5,11:16))];
DATA_TC = [cell2mat(Cell_TC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_TC(3:5,11:16))];
DATA_DMSO = [cell2mat(DMSO_TC(9:11,2:9)),[nan;nan;nan],cell2mat(DMSO_TC(9:11,11:16))];
DATA_DMSO_avg = mean(DATA_DMSO,1);

DATA    = zeros(NT,NC,NR);
DATA_sc = zeros(NT,NC,NR);
for i = 1:NT
    k = i;
    DATA(i,:,:) = DATA_DMSO_avg(3*k-2:3*k).*DATA_TC(:,3*k-2:3*k)/100;
    DATA_sc(i,:,:) = squeeze(DATA(i,:,:)).*DATA_SC(:,3*k-2:3*k)/100;
end
% DATA(1,:,:) = init;
% DATA_sc(1,:,:) = 0.25/100*init;



r = squeeze(DATA_sc(1,:,:))./squeeze(DATA(1,:,:));

ratio = mean(r,2,'omitnan');

DATA(:,:,3) = [];
DATA_sc(:,:,3) = [];
NR = 2;




%% Optimization bound (Hill2_switching_death)
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
E_ub  = 100;
m_lb  = 0;
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
beq = [0.25/100,1-0.25/100];
% Aeq = [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
% beq = [1];

num_optim= 100;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(2) = xi(3) + rand*lam_s_range(2);
    xi(4) = rand*lam_d_range(2);
    xi(12) = xi(13);
    xi(1)  = 0.25/100;
    xi(11) = 1-0.25/100;
    x_init = [x_init;xi];
end





%% Optimization (Point estimate death)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_alt_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd) + get_like_onepop_h(DATA_sc,x,Time,Conc,NR,NC,NT,s,cmd,1);
fval_hist_pe_DIP   = [];
params_hist_pe_DIP = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,[],[],lb,ub,[],options1); 
    fval_hist_pe_DIP = [fval_hist_pe_DIP,ff];
    params_hist_pe_DIP = [params_hist_pe_DIP;xx];
%     catch
%     end
end
t = toc
[of_DIP,oi] = min(fval_hist_pe_DIP);
opt_xx_pe_DIP = params_hist_pe_DIP(oi,:);


%% nDIP
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
E_ub  = 100;
m_lb  = 0;
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
beq = [0.25/100,1-0.25/100];
% Aeq = [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
% beq = [1];

num_optim= 100;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(2) = xi(3) + rand*lam_s_range(2);
    xi(4) = rand*lam_d_range(2);
    xi(12) = xi(13);
    xi(1)  = 0.25/100;
    xi(11) = 1-0.25/100;
    x_init = [x_init;xi];
end






%% Optimization (Point estimate death)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_alt_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd) + get_like_onepop_h(DATA_sc,x,Time,Conc,NR,NC,NT,s,cmd,1);
fval_hist_pe_nDIP   = [];
params_hist_pe_nDIP = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,[],[],lb,ub,[],options1); 
    fval_hist_pe_nDIP = [fval_hist_pe_nDIP,ff];
    params_hist_pe_nDIP = [params_hist_pe_nDIP;xx];
%     catch
%     end
end
t = toc
[of_nDIP,oi] = min(fval_hist_pe_nDIP);
opt_xx_pe_nDIP = params_hist_pe_nDIP(oi,:);




%% Compute the AIC


nDIP_AIC = 2*10 + 2*of_nDIP;
DIP_AIC = 2*12 + 2*of_DIP;

    

save_name = strcat('Result/In_vitro_Time_all_h_TCSC_AIC.mat');

save(save_name)


