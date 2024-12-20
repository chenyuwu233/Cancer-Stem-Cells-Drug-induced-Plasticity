parpool('local',20)
warning('off','MATLAB:integral:NonFiniteValue')



%%  Setting the Concentration and Time points

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


DATA    = zeros(NT,NC,NR);

DATA(3,:,1) = DATA_TC(:,1)*mean(DMSO_1)/100;
DATA(3,:,2) = DATA_TC(:,2)*mean(DMSO_2)/100;


load('Forward_D24_v.mat')
DATA(2,:,1) = DATA_rate(1,:);
DATA(2,:,2) = DATA_rate(2,:);

load('T12_2T_v.mat')
DATA(1,:,1) = T12_2T(1);
DATA(1,:,2) = T12_2T(2);

DATA = round(DATA);





%% Optimization bound (DIP) TD2



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
func = @(x) get_like_TD_meanflow_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_TD2   = [];
params_hist_pe_TD2 = [];
tic
parfor j = 1:num_optim
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_TD2 = [fval_hist_pe_TD2,ff];
    params_hist_pe_TD2 = [params_hist_pe_TD2;xx];
end
t_TD2 = toc
[of_TD2,oi] = min(fval_hist_pe_TD2);
opt_xx_pe_TD2 = params_hist_pe_TD2(oi,:);





%% Logistic time delay with 2 parameters 







    
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
func = @(x) get_like_TD_meanflow_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_nTD2   = [];
params_hist_pe_nTD2 = [];
tic
parfor j = 1:num_optim
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_nTD2 = [fval_hist_pe_nTD2,ff];
    params_hist_pe_nTD2 = [params_hist_pe_nTD2;xx];
end
t_nTD2 = toc
[of_nTD2,oi] = min(fval_hist_pe_nTD2);
opt_xx_pe_nTD2 = params_hist_pe_nTD2(oi,:);


%% Compute the AIC



nDIP_TD2_num_par = 8;
DIP_TD2_num_par = 14;


nDIP_TD2_AIC = 2*nDIP_TD2_num_par + 2*of_nTD2;
DIP_TD2_AIC = 2*DIP_TD2_num_par + 2*of_TD2;



nDIP_TD2_AICc = nDIP_TD2_AIC + (2*nDIP_TD2_num_par^2 + 2*nDIP_TD2_num_par)/(3*9*2 - nDIP_TD2_num_par - 1);
DIP_TD2_AICc = DIP_TD2_AIC + (2*DIP_TD2_num_par^2 + 2*DIP_TD2_num_par)/(3*9*2 - DIP_TD2_num_par - 1);




    

save_name = strcat('Result/In_vitro_Gastric_CNSC_CSC_TD.mat');

save(save_name)


