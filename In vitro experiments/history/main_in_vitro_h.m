%%  Generating the theta:
% - s: sub-population
% - p: Initial proportion
% - 
% % 


parpool('local',20)
warning('off','MATLAB:integral:NonFiniteValue')



%%  Setting the Concentration and Time points

% init   = 1e5;
s      = 2;
Conc = [0,4,8];
NT   = 4;
Time = [0,4,22,46];
NC   = length(Conc);
NR   = 3;
cmd  = 'CSC_DIS';
lam_s_range       = [0,0.1];
lam_d_range       = [0,0.5];



%% Load data

Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figures 2A');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 2B');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 2B');


DATA_SC = [cell2mat(Cell_SC(3:5,2:7)),cell2mat(Cell_SC(3:5,12:17))];
DATA_TC = [cell2mat(Cell_TC(3:5,2:7)),cell2mat(Cell_TC(3:5,11:16))];
DATA_DMSO = [cell2mat(DMSO_TC(9:11,2:7)),cell2mat(DMSO_TC(9:11,11:16))];
DATA_DMSO_avg = mean(DATA_DMSO,1);

DATA    = zeros(NT,NC,NR);
for i = 1:4
    DATA(i,:,:) = DATA_DMSO_avg(3*i-2:3*i).*DATA_TC(:,3*i-2:3*i)/100;
end







%% Optimization bound (Hill2_switching_death)
%  theta = [{alpha,beta,nu,b_beta,E_beta}_s,c]


alpha_lb = 0;
alpha_ub = 10-1e-6;
alpha_d_ub = 0;
beta_lb  = 1e-6;
beta_ub  = 10;
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
     0,0,0,1,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0];
b = [lam_s_range(2);lam_d_range(2)];

Aeq = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
beq = [0.25/100,1-0.25/100];

num_optim= 20;


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





    




% %% Optimization parts (CI)
% 
% options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
% 
% 
% 
% 
% 
% 
% 
% 
% B_num = 13;
% B_sample = 100;
% B_parameter = [];
% 
% 
% parfor i = 1:B_sample
% 
%     Data_i =zeros(NT,NC,B_num);
%     for j = 1:B_num
%         ri = randi(NR);
%         Data_i(:,:,j) =  DATA(:,:,ri);
%     end
%     
%     func = @(x) get_like(Data_i,x,Time,Conc,B_num,NC,NT,s,cmd);
%     
%     fval_hist   = [];
%     params_hist = [];
%     
%     for j = 1:num_optim
%         [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),A,b,[],[],lb,ub,[],options1); 
%         fval_hist = [fval_hist,ff];
% %         fval_hist = [fval_hist,func(xx)];
%         params_hist = [params_hist;xx];
%     end
%     [of,oi] = min(fval_hist);
%     opt_xx = params_hist(oi,:);
%     B_parameter = [B_parameter;opt_xx];
%     fprintf('Loop %d complete\n',i)
% end


%% Optimization (Point estimate death)
options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
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



    

save_name = strcat('Result/In_vitro_h.mat');

save(save_name)


