% Code for generating Bootstrap confidence interval estimation in base
% experiment. Note that this script requires long time to run.

%%  Generating the theta:
parpool('local',100)
warning('off','MATLAB:integral:NonFiniteValue')


%%

seed_num = 130;

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
nu_diff_range = [0.05,0.4];
b_beta_range    = [0.8,0.9];
b_nu_range      = [1+1e-6,1.1];
E_range         = [0.0625,2.5];
c_range         = [0,10];


c      = rand*10;
% p      = rand*(p_range(2)-p_range(1)) + p_range(1);
% alpha1 = rand*(alpha_range(2)-alpha_range(1)) + alpha_range(1);
% alpha2 = rand*(alpha_range(2)-alpha_range(1)) + alpha_range(1);
beta1  = rand*(beta_s_range(2)-beta_s_range(1)) + beta_s_range(1);
beta2  = rand*(beta_d_range(2)-beta_d_range(1)) + beta_d_range(1);
alpha1 = beta1 + rand*lam_s_range(2);
alpha2 = beta2 + rand*lam_d_range(2);
nu12   = rand*nu_diff_range(2);
nu21   = 0;
b1_beta = 1;
b2_beta = rand*(b_beta_range(2)-b_beta_range(1)) + b_beta_range(1);
b2_nu   = rand*(b_nu_range(2)-b_nu_range(1)) + b_nu_range(1);
E2_beta = rand*(E_range(2)-E_range(1)) + E_range(1);
E2_nu = rand*(E_range(2)-E_range(1)) + E_range(1);

theta   = [alpha1,beta1,nu12,1,1,1,1,...
              alpha2,beta2,nu21,b2_beta,E2_beta,b2_nu,E2_nu,c]

A = [alpha1-beta1,nu12;0,alpha2-beta2];


get_stable_p(A)




%%


tic
DATA = Switching_gen(init,theta,Time,Conc,NR,NC,NT,s,cmd);
t_gen = toc

tic
opt_fval = get_like(DATA,theta,Time,Conc,NR,NC,NT,s,cmd);
t_like = toc






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
b_beta_lb  = 0.5;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 1.5;
E_lb  = 1e-6;
E_ub  = 5;
c_lb  = 0;
c_ub  = 10;

lb = [alpha_lb,beta_lb,nu_lb,1,1,1,1,...
        alpha_lb,beta_lb,nu_lb,b_beta_lb,E_lb,b_nu_lb,E_lb,c_lb];



ub = [alpha_ub,beta_ub,nu_ub,1,1,1,1,...
        alpha_ub,beta_ub,nu_d_ub,b_beta_ub,E_ub,b_nu_ub,E_ub,c_ub];


A = [1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0;
     -1,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0];
b = [lam_s_range(2);lam_d_range(2);0;0];

num_optim= 20;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(1) = xi(2) + rand*lam_s_range(2);
    % xi(8) = xi(9);
    xi(3) = rand*lam_d_range(2);
    x_init = [x_init;xi];
end






    




%% Optimization parts (CI)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');








B_num = 13;
B_sample = 100;
B_parameter = [];


parfor i = 1:B_sample

    Data_i =zeros(NT,NC,B_num);
    for j = 1:B_num
        ri = randi(NR);
        Data_i(:,:,j) =  DATA(:,:,ri);
    end

    func = @(x) get_like(Data_i,x,Time,Conc,B_num,NC,NT,s,cmd);

    fval_hist   = [];
    params_hist = [];

    for j = 1:num_optim
        [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,[],[],lb,ub,[],options1); 
        fval_hist = [fval_hist,ff];
%         fval_hist = [fval_hist,func(xx)];
        params_hist = [params_hist;xx];
    end
    [of,oi] = min(fval_hist);
    opt_xx = params_hist(oi,:);
    B_parameter = [B_parameter;opt_xx];
    fprintf('Loop %d complete\n',i)
end








    

save_name = strcat('Result_CI/Boot_CI_CSC_DIS_alt_',num2str(seed_num),'.mat');

save(save_name)


