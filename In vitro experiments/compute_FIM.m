


load('Result/In_vitro_AIC_BaF3_21_YY_ip.mat')

% parpool('local',40)
% warning('off','MATLAB:integral:NonFiniteValue')


%%

options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','iter','algorithm','sqp');
[xx,ff,~,out,~,g,H_opt]  = fmincon(func,opt_xx_pe_DIP_Asy,A,b,Aeq,beq,lb,ub,[],options1); 


%%



theta = opt_xx_pe_DIP_Asy;

N = 400;
init = round(mean(DATA(1,:,:),"all",'omitnan'));

H_N = zeros(length(theta),length(theta));

parfor num = 1:N

    %% Gen data
    
    
    
    tic
    DATA_n = Switching_gen_h(init, theta, Time, Conc, NR, NC, NT, s,cmd);
    td = toc
    
    %% Compute M hessian estimates
    M = 10;
    
    H_n = zeros([length(theta),length(theta)]);
    
    func_n = @(x) get_like_BaF3(DATA_n,x,Time,Conc,NR,NC,NT,s,cmd);
    tic
    for m = 1:M
        
        
        Delta = zeros(size(theta));
        c     = 1e-6;
        
        for i = 1:length(Delta)
            if rand < 0.5
                Delta(i) = c;
            else
                Delta(i) = -c;
            end
        end
    
        theta_PD = theta + Delta;
        theta_MD = theta - Delta;
    
        Eps = zeros(size(theta));
        c     = 1e-6;
        
        for i = 1:length(Eps)
            if rand < 0.5
                Eps(i) = c;
            else
                Eps(i) = -c;
            end
        end
    
        theta_PD_PE = theta_PD + Eps;
        theta_PD_ME = theta_PD - Eps;
        theta_MD_PE = theta_MD + Eps;
        theta_MD_ME = theta_MD - Eps;
    
        G_P = (func_n(theta_PD_PE) - func_n(theta_PD_ME))/2./Eps;
        G_M = (func_n(theta_MD_PE) - func_n(theta_MD_ME))/2./Eps;
        delta_G = G_P - G_M;
        G_H = delta_G'*(Delta.^(-1))./2;
        H_m = 1/2*(G_H + G_H');
        
        H_n = H_n + H_m;
        
    end
    tavg = toc

    H_N = H_N + H_n./M;

end

H_MC = H_N./N;

FIM_N = H_MC;



save('Result/In_vitro_AIC_BaF3_21_YY_ip_FIM.mat')





