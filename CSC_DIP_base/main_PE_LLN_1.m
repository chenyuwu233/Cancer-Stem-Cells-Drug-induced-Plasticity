%%  Generating the theta:
% parpool('local',20)
% warning('off','MATLAB:integral:NonFiniteValue')


%%  Access the existing data



for ii = 31:40
seed_num = ii;

rng(seed_num)


data_name = strcat('Result_GS/CIPE_CSC_DIS_',num2str(seed_num),'.mat');
load(data_name)




ub(10) = 0;
ub(end) = 1e4;


    


%% Optimization (Point estimate)

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off','algorithm','sqp');
func = @(x) get_like_LLM(DATA,x,Time,Conc,NR,NC,NT,s,cmd);
fval_hist_pe   = [];
params_hist_pe = [];
H_history = cell(1,num_optim);
tic
parfor j = 1:num_optim
    try
    [xx,ff,~,out,~,g,H]  = fmincon(func,x_init(j,:),A,b,[],[],lb,ub,[],options1); 
    fval_hist_pe = [fval_hist_pe,ff];
    params_hist_pe = [params_hist_pe;xx];
    H_history{j} = H;
    catch
    end
    j
end
t = toc
[of,oi] = min(fval_hist_pe);
opt_xx_pe = params_hist_pe(oi,:);
opt_Hessian = H_history{oi};


%% Fmincon check point
wrong_idx = 0;
fmin_diff_hist = [];
for i = 1:size(params_hist_pe,1)
    fmin_diff_hist = [fmin_diff_hist,abs(func(params_hist_pe(i,:)) - fval_hist_pe(i))];
    if abs(func(params_hist_pe(i,:)) - fval_hist_pe(i)) >= 1e-3
        wrong_idx = 1;
    end
end



    

save_name = strcat('Result_GS/LLM_Results/PE_CSC_DIS_',num2str(seed_num),'.mat');

save(save_name)

end



%% 

for ii = 31:130
    seed_num = ii;
    data_name = strcat('Result_GS/LLM_Results/PE_CSC_DIS_',num2str(seed_num),'.mat');
    load(data_name)
    theta_lln = theta;
    lb = theta_lln;
    ub = theta_lln;
    lb(end) = 0;
    ub(end) = 1e4;
    [xx,ff] = fmincon(func,theta_lln,[],[],[],[],lb,ub,[],options1);
    if ff < of
        keyboard
    end
end