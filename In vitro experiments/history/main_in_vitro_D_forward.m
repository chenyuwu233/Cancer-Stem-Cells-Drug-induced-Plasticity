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
NT   = 2;
Time = [0,36];
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

% load('T12_2T.mat')
% DATA(1,:,1) = T12_2T;
% DATA(1,:,2) = T12_2T;

% load('rate12.mat')
% DATA(1,:,1) = mean(DMSO_1)/rate_12;
% DATA(1,:,2) = mean(DMSO_2)/rate_12;

DATA(2,:,1) = DATA_TC(:,1)*mean(DMSO_1)/100;
DATA(2,:,2) = DATA_TC(:,2)*mean(DMSO_2)/100;


% load('rate24.mat')
% DATA(1,:,1) = DATA(2,:,1)./rate;
% DATA(1,:,2) = DATA(2,:,2)./rate;

load('Forward_D24.mat')
DATA(1,:,1) = DATA_rate;
DATA(1,:,2) = DATA_rate;




%% Optimization bound (Hill2_switching_death)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]


% alpha_lb = 0;
% alpha_ub = 1-1e-6;
% alpha_d_ub = 0;
% beta_lb  = 1e-6;
% beta_ub  = 1;
% nu_lb = 0;
% nu_ub = 1-1e-6;
% nu_d_ub = 0;
% b_beta_lb  = 0.5;
% b_beta_ub  = 1;
% b_nu_lb    = 1;
% b_nu_ub    = 1.5;
% E_lb  = 1e-6;
% E_ub  = 16;
% c_lb  = 0;
% c_ub  = 300;
% 
% lb = [0,alpha_lb,beta_lb,nu_lb,1,1,1,1,...
%         0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,b_nu_lb,E_lb,c_lb];
% 
% 
% 
% ub = [1,alpha_ub,beta_ub,nu_ub,1,1,1,1,...
%         1,alpha_ub,beta_ub,0,b_beta_ub,E_ub,b_nu_ub,E_ub,c_ub];
% 
% 
% A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%      0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% b = [0.1;0];
% 
% Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%        0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0];
% beq = [0,1];
% 
% num_optim= 20;
% 
% 
% x_init = [];
% for i = 1:num_optim
%     xi = rand(1,length(ub)).*(ub-lb)+lb;
%     xi(2) = xi(3) + rand*lam_s_range(2);
%     xi(4) = rand*lam_d_range(2);
%     xi(10) = xi(11);
%     x_init = [x_init;xi];
% end






%% Optimization bound (Hill2_switching_death)
%  theta = [{alpha,beta,nu,b_beta,E_beta,n_beta}_s,c]


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
m_lb  = 0;
m_ub  = 5;
c_lb  = 0;
c_ub  = 1000;

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
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd) + (x(14) - 10.2182)^2 ;
% func = @(x) get_like_alt(DATA,x,Time,Conc,NR,NC,NT,s,cmd) + 1/0.0286^2*(x(2)-x(3)-0.0723)^2;
% func = @(x) get_like_alt_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd) + (x(16) - 8.4676)^2 + (x(17) - 1.8479)^2 + (x(2) - x(4))^2 + 1.6*(x(15) - 0.9564)^2;
func = @(x) get_like_alt_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
%     catch
%     end
end
t = toc
[of,oi] = min(fval_hist_pe);
opt_xx_pe = params_hist_pe(oi,:);


%% Fmincon check point
wrong_idx = 0;
fmin_diff_hist = [];
for i = 1:size(params_hist_pe,1)
    fmin_diff_hist = [fmin_diff_hist,abs(func(params_hist_pe(i,:)) - fval_hist_pe(i))];
    if abs(func(params_hist_pe(i,:)) - fval_hist_pe(i)) >= 1e-3
        wrong_idx = 1;
    end
end



    

save_name = strcat('Result/In_vitro_Conc_2448_h_forward.mat');

save(save_name)


