%%  Generating the theta:
% - s: sub-population
% - p: Initial proportion
% - 
% % 


parpool('local',50)
warning('off','MATLAB:integral:NonFiniteValue')


%% Load data


load("Mixture_Data.mat")


DATA = permute(BF_21,[3,2,1]); % permute the data matrix to [Time, Conc, replicates]




%%  Setting the Concentration and Time points

% init   = 5000;
s      = 2;
[NT,NC,NR] = size(DATA);
Time = [0 : NT-1]*3;
Conc = 10^(6)*[0  31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6) ];

cmd  = 'CSC_DIS';








%% Optimization (nDIP_nAsy)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]


alpha_lb = 0;
alpha_ub = 1-1e-6;
beta_lb  = 1e-6;
beta_ub  = 1;
b_beta_lb  = 0.878;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 50;
m_lb  = 1e-3;
m_ub  = 20;
c_lb  = 0;
c_ub  = 100;

lb = [0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,1,1,1,...
      0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,1,1,1,c_lb];



ub = [1,alpha_ub,beta_ub,0,b_beta_ub,E_ub,m_ub,1,1,1,...
      1,alpha_ub,beta_ub,0,b_beta_ub,E_ub,m_ub,1,1,1,c_ub];




A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0];
b = [0.06;0;0.06;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = [2/3;1/3];
% Aeq = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
% beq = [1];

num_param_nDIP_nAsy = sum(lb ~= ub) - length(beq);

num_optim= 25;


x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(11) = 1 - xi(1);
    xi(2) = xi(3) + rand*0.06;
    xi(12) = xi(13) + rand*0.06;
    x_init = [x_init;xi];
end
    


options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_BaF3(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_nDIP_nAsy   = [];
params_hist_pe_nDIP_nAsy = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_nDIP_nAsy = [fval_hist_pe_nDIP_nAsy,ff];
    params_hist_pe_nDIP_nAsy = [params_hist_pe_nDIP_nAsy;xx];
%     catch
%     end
end
t_nDIP_nAsy = toc
[of_nDIP_nAsy,oi] = min(fval_hist_pe_nDIP_nAsy);
opt_xx_pe_nDIP_nAsy = params_hist_pe_nDIP_nAsy(oi,:);


%% Optimization (nDIP_Asy)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]


alpha_lb = 0;
alpha_ub = 1-1e-6;
beta_lb  = 1e-6;
beta_ub  = 1;
b_beta_lb  = 0.878;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 50;
m_lb  = 1e-3;
m_ub  = 20;
c_lb  = 0;
c_ub  = 100;

lb = [0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,1,1,1,...
      0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,1,1,1,c_lb];



ub = [1,alpha_ub,beta_ub,0.1,b_beta_ub,E_ub,m_ub,1,1,1,...
      1,alpha_ub,beta_ub,0.1,b_beta_ub,E_ub,m_ub,1,1,1,c_ub];


A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0];
b = [0.06;0;0.06;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = [2/3;1/3];
% Aeq = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
% beq = [1];

num_param_nDIP_Asy = sum(lb ~= ub) - length(beq);



x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(11) = 1 - xi(1);
    xi(2) = xi(3) + rand*0.06;
    xi(12) = xi(13) + rand*0.06;
    x_init = [x_init;xi];
end
    


options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','iter','algorithm','sqp');
func = @(x) get_like_BaF3(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_nDIP_Asy   = [];
params_hist_pe_nDIP_Asy = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_nDIP_Asy = [fval_hist_pe_nDIP_Asy,ff];
    params_hist_pe_nDIP_Asy = [params_hist_pe_nDIP_Asy;xx];
%     catch
%     end
end
t_nDIP_Asy = toc
[of_nDIP_Asy,oi] = min(fval_hist_pe_nDIP_Asy);
opt_xx_pe_nDIP_Asy = params_hist_pe_nDIP_Asy(oi,:);



%% Optimization (DIP_nAsy)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]


alpha_lb = 0;
alpha_ub = 1-1e-6;
beta_lb  = 1e-6;
beta_ub  = 1;
b_beta_lb  = 0.878;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 50;
m_lb  = 1e-3;
m_ub  = 20;
c_lb  = 0;
c_ub  = 100;

lb = [0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,b_nu_lb,E_lb,m_lb,...
      0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,b_nu_lb,E_lb,m_lb,c_lb];



ub = [1,alpha_ub,beta_ub,0,b_beta_ub,E_ub,m_ub,b_nu_ub,E_ub,m_ub,...
      1,alpha_ub,beta_ub,0,b_beta_ub,E_ub,m_ub,b_nu_ub,E_ub,m_ub,c_ub];


A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0];
b = [0.06;0;0.06;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = [2/3;1/3];

% Aeq = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
% beq = [1];

num_param_DIP_nAsy = sum(lb ~= ub) - length(beq);




x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(11) = 1 - xi(1);
    xi(2) = xi(3) + rand*0.06;
    xi(12) = xi(13) + rand*0.06;
    x_init = [x_init;xi];
end
    


options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_BaF3(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_DIP_nAsy   = [];
params_hist_pe_DIP_nAsy = [];
tic
parfor j = 1:num_optim
%     try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_DIP_nAsy = [fval_hist_pe_DIP_nAsy,ff];
    params_hist_pe_DIP_nAsy = [params_hist_pe_DIP_nAsy;xx];
%     catch
%     end
end
t_DIP_nAsy = toc
[of_DIP_nAsy,oi] = min(fval_hist_pe_DIP_nAsy);
opt_xx_pe_DIP_nAsy = params_hist_pe_DIP_nAsy(oi,:);


Optimization (DIP_Asy)
 theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]


alpha_lb = 0;
alpha_ub = 1-1e-6;
beta_lb  = 1e-6;
beta_ub  = 1;
b_beta_lb  = 0.878;
b_beta_ub  = 1;
b_nu_lb    = 1;
b_nu_ub    = 2;
E_lb  = 1e-6;
E_ub  = 50;
m_lb  = 1e-3;
m_ub  = 20;
c_lb  = 0;
c_ub  = 100;

lb = [0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,b_nu_lb,E_lb,m_lb,...
      0,alpha_lb,beta_lb,0,b_beta_lb,E_lb,m_lb,b_nu_lb,E_lb,m_lb,c_lb];



ub = [1,alpha_ub,beta_ub,0.1,b_beta_ub,E_ub,m_ub,b_nu_ub,E_ub,m_ub,...
      1,alpha_ub,beta_ub,0.1,b_beta_ub,E_ub,m_ub,b_nu_ub,E_ub,m_ub,c_ub];


A = [0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0];
b = [0.06;0;0.06;0];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = [2/3;1/3];
% Aeq = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
% beq = 1;

num_param_DIP_Asy = sum(lb ~= ub) - length(beq);




x_init = [];
for i = 1:num_optim
    xi = rand(1,length(ub)).*(ub-lb)+lb;
    xi(11) = 1 - xi(1);
    xi(2) = xi(3) + rand*0.06;
    xi(12) = xi(13) + rand*0.06;
    x_init = [x_init;xi];
end



options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_BaF3(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe_DIP_Asy   = [];
params_hist_pe_DIP_Asy = [];
tic
parfor j = 1:num_optim
    try
    [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,Aeq,beq,lb,ub,[],options1); 
    fval_hist_pe_DIP_Asy = [fval_hist_pe_DIP_Asy,ff];
    params_hist_pe_DIP_Asy = [params_hist_pe_DIP_Asy;xx];
    catch
    end
end
t_DIP_Asy = toc
[of_DIP_Asy,oi] = min(fval_hist_pe_DIP_Asy);
opt_xx_pe_DIP_Asy = params_hist_pe_DIP_Asy(oi,:);

%% Compute the AIC


nDIP_nAsy_AIC = 2*num_param_nDIP_nAsy + 2*of_nDIP_nAsy;
nDIP_Asy_AIC = 2*num_param_nDIP_Asy + 2*of_nDIP_Asy;
DIP_nAsy_AIC = 2*num_param_DIP_nAsy + 2*of_DIP_nAsy;
DIP_Asy_AIC = 2*num_param_DIP_Asy + 2*of_DIP_Asy;




save_name = strcat('Result/In_vitro_AIC_BaF3_21_YY_ip.mat');

save(save_name)


