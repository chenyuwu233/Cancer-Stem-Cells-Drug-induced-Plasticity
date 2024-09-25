%%  Generating the theta:
% - s: sub-population
% - p: Initial proportion
% - 
% 
parpool('local',20)
warning('off','MATLAB:integral:NonFiniteValue')
% addpath('C:\Users\euclid\OneDrive\UMN-My-gear\Phenotypic_Switching_model\Switching_model(MATLAB)\Main functions')

seed_num = 72;

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
lam_d_range       = [0,0.5];
nu_dediff_range = [0,0.05];
b_beta_range    = [0.8,0.9];
b_nu_range      = [1+1e-6,1.1];
E_range         = [0.0625,2.5];
c_range         = [0,10];


c      = rand*10;
% p      = rand;
% alpha1 = rand*(alpha_range(2)-alpha_range(1)) + alpha_range(1);
% alpha2 = rand*(alpha_range(2)-alpha_range(1)) + alpha_range(1);
beta1  = rand*(beta_s_range(2)-beta_s_range(1)) + beta_s_range(1);
beta2  = rand*(beta_d_range(2)-beta_d_range(1)) + beta_d_range(1);
alpha1 = beta1 + rand*lam_s_range(2);
alpha2 = beta2;
nu12   = rand*lam_d_range(2);
nu21   = 0;
b1_beta = 1;
b2_beta = rand*(b_beta_range(2)-b_beta_range(1)) + b_beta_range(1);
b2_nu   = rand*(b_nu_range(2)-b_nu_range(1)) + b_nu_range(1);
E2_beta = rand*(E_range(2)-E_range(1)) + E_range(1);
E2_nu = rand*(E_range(2)-E_range(1)) + E_range(1);

theta   = [alpha1,beta1,nu12,1,1,1,1,...
              alpha2,beta2,nu21,b2_beta,E2_beta,b2_nu,E2_nu,c]
theta_null = [alpha1,beta1,nu12,1,1,1,1,...
              0,beta2,nu21,b2_beta,E2_beta,b2_nu,E2_nu,c]

A = [alpha1-beta1,nu12;0,-beta2];


eig(A)

get_stable_p(A)



%% cap test

% no_cap_cap_hist = [];
% cap_cap_hist = [];
% cap_null_hist = [];
% null_null_hist = [];
% 
% for cap = 0:6
% 
%     tic
%     DATA_no_cap = Switching_gen(init,theta,Time,Conc,NR,NC,NT,s,cmd);
%     t_gen1 = toc
%     
%     tic
%     DATA_cap = Switching_gen(init,theta,Time,Conc,NR,NC,NT,s,cmd,cap);
%     t_gen2 = toc
%     
%     tic
%     DATA_null = Switching_gen(init,theta_null,Time,Conc,NR,NC,NT,s,cmd);
%     t_gen3 = toc
%     
%     tic
%     opt_fval_no_cap_cap = get_like(DATA_no_cap,theta,Time,Conc,NR,NC,NT,s,cmd)
%     no_cap_cap_hist = [no_cap_cap_hist,opt_fval_no_cap_cap];
%     t_like = toc
%     
%     tic
%     opt_fval_cap_cap = get_like(DATA_cap,theta,Time,Conc,NR,NC,NT,s,cmd)
%     cap_cap_hist = [cap_cap_hist,opt_fval_cap_cap];
%     t_like = toc
%     
%     tic
%     opt_fval_cap_null = get_like(DATA_cap,theta_null,Time,Conc,NR,NC,NT,s,cmd)
%     cap_null_hist = [cap_null_hist,opt_fval_cap_null];
%     t_like = toc
%     
%     tic
%     opt_fval_null_null = get_like(DATA_null,theta_null,Time,Conc,NR,NC,NT,s,cmd)
%     null_null_hist = [null_null_hist,opt_fval_null_null];
%     t_like = toc
% 
% end




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
     0,0,1,0,0,0,0,1,-1,0,0,0,0,0,0];
b = [lam_s_range(2);lam_d_range(2)];

num_optim= 20;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(1) = xi(2) + rand*lam_s_range(2);
    xi(8) = xi(9);
    xi(3) = rand*lam_d_range(2);
    x_init = [x_init;xi];
end

%%

opt_xx_pe_hist = [];
for cap = 0:3
    


    tic
    DATA = Switching_gen(init,theta,Time,Conc,NR,NC,NT,s,cmd,cap);
    t_gen = toc

    tic
    opt_fval_cap_null = get_like(DATA,theta_null,Time,Conc,NR,NC,NT,s,cmd);
    t_like_null = toc

    tic
    opt_fval_cap = get_like(DATA,theta,Time,Conc,NR,NC,NT,s,cmd);
    t_like = toc
    
%% 

    tic
    DATA_null = Switching_gen(init,theta_null,Time,Conc,NR,NC,NT,s,cmd);
    t_gen = toc
    
    
 
    
    
    %% Optimization (Point estimate)
    
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
    func = @(x) get_like(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
    fval_hist_pe   = [];
    params_hist_pe = [];
    
    parfor j = 1:num_optim
        try
        [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,[],[],lb,ub,[],options1); 
        fval_hist_pe = [fval_hist_pe,ff];
        params_hist_pe = [params_hist_pe;xx];
        catch
        end
    end
    [of,oi] = min(fval_hist_pe);
    opt_xx_pe = params_hist_pe(oi,:);
    opt_xx_pe_hist = [opt_xx_pe_hist;opt_xx_pe];


end


% %% Fmincon check point
% wrong_idx = 0;
% fmin_diff_hist = [];
% for i = 1:size(params_hist_pe,1)
%     fmin_diff_hist = [fmin_diff_hist,abs(func(params_hist_pe(i,:)) - fval_hist_pe(i))];
%     if abs(func(params_hist_pe(i,:)) - fval_hist_pe(i)) >= 1e-3
%         wrong_idx = 1;
%     end
% end

%%

    

save_name = strcat('Result/PE_CSC_DIS_a12_',num2str(seed_num),'.mat');

save(save_name)


