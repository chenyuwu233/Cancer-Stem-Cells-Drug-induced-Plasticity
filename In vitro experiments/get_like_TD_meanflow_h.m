%% Description:
%  This is a function to obtain the Covaraince matrix for swithing model
%  Important assumption: the time interval is identical, theta contains the
%  initial proportion.
%
%  Required function:
%  get_path_like.m
%  get_Cov_mat.m
%  get_Cov.m
%  get_sig_j.m
%  get_m.m
%  
%  Input:
%  - DATA: NT X NC X NR tensor that records all the data (In this situation, we assume data have unify initial cell count)
%  - theta: 1 x s*(s+4)+1 that includes the initial proportion, birth rate, death rate, and
%           switching rate among every sub-types, the last element is drug effect.
%                   [{p,alhpa,beta,{nu}_{s-1},b,E,(n)}_s,c]
%  - Time: 1 x NT vector that records all the time points, Time points does
%  not need to be equal-distance.
%  - Conc: 1,2 x NC vector that records all the concentration points
%  - NR: number of replicates
%  - NC: number of concentration points
%  - NT: number of time points
%  - s:  number of distinct sub-types
%  - cmd:   string command about what model we used for drug dependency:
%      - 'Linear_switching': Assume the switching rate (nu)_i 'Linearly'
%      depends on the drug concentration with parameters (b)_i
%
%      - 'Single_Broad': Assume the switching rate (nu)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i
%
%      - 'Single_Broad_d': Assume the death rate (beta)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i
%
%      - 'Hill2_swtiching_2': Assume the switching rate (nu)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently.
%
%      - 'Hill2_death_2': Assume the death rate (beta)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently.
%
%      - 'Multi_Selective': Assume the switching rate (nu)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently. We applied the drug independently.
%
%      - 'Multi_Selective_d': Assume the death rate (beta)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently. We applied the drug independently.
%
%      - 'Hill2_switching_death': Assume the switching rate (nu)_i and
%      death rate (beta)_i 'Hilly' depends on the drug concentration with
%      parameters (b_nu,E_nu,b_beta,E_beta)_i. 
%      Note that this model only has one type of drug.
%
%      - 'Hill2_switching_Hill2_death': Assume the switching rate (nu)_i
%      and death rate (beta)_i 'Hilly' depends on the drug concentration
%      with parameters (b_nu,E_nu,b_beta,E_beta)_i. 
%      Note that this model allows two type of drug, one affect the
%      switching rate (nu)_i and one affect the death rate (beta)_i.
%      
%      - 'Hill2_switching_Hill2_switching': Assume the switching rate (nu)_i
%      'Hilly' depends on the drug concentration
%      with parameters (b_nu1,E_nu1,b_nu2,E_nu2)_i. 
%      Note that this model allows two type of drug, one affect the
%      switching rate with (b_nu1)_i and one affect the switching rate with
%      (b_nu2)_i
%
%      - 'Hill2_death_Hill2_death': Assume the death rate (beta)_i 
%      'Hilly' depends on the drug concentration
%      with parameters (b_beta1,E_beta1,b_beta2,E_beta2)_i. 
%      Note that this model allows two type of drug, one affect the
%      death rate with (b_beta1)_i and one affect the death rate with
%      (b_beta2)_i
%  
%      - 'Hill2_switching_Hill2_death_i': Assume the switching rate (nu)_i
%      and death rate (beta)_i 'Hilly' depends on the drug concentration
%      with parameters (b_nu,E_nu,b_beta,E_beta)_i. 
%      Note that this model allows two type of drug, one affect the
%      switching rate (nu)_i and one affect the death rate (beta)_i. We
%      applied the drug independently.
%      
%      - 'Multi_Broad': Assume the switching rate (nu)_i
%      'Hilly' depends on the drug concentration
%      with parameters (b_nu1,E_nu1,b_nu2,E_nu2)_i. 
%      Note that this model allows two type of drug, one affect the
%      switching rate with (b_nu1)_i and one affect the switching rate with
%      (b_nu2)_i. We applied the drug independently.
%
%      - 'Hill2_death_Hill2_death_i': Assume the death rate (beta)_i 
%      'Hilly' depends on the drug concentration
%      with parameters (b_beta1,E_beta1,b_beta2,E_beta2)_i. 
%      Note that this model allows two type of drug, one affect the
%      death rate with (b_beta1)_i and one affect the death rate with
%      (b_beta2)_i. We applied the drug independently.
%
%      - 'CSC_DIS': One drug which can increase death rate and switching rate of the DC.
%        In this experiment we have the following parameters set:
%        θ=[α_1,β_1,ν_12, α_2=0,β_2,ν_21,b_beta_2,E_beta_2,b_nu_2,E_nu_2,tv,c]
%
%  Output:
%  The Likelihood of the DATA


function ret = get_like_TD_meanflow_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd)
    theta = x;
    c = theta(end)';
    

    theta(end) = [];
    tv = theta(end-1:end); % Variable to control the delay effect 
    theta(end-1:end) = [];


    Theta = reshape(theta,[],s)';
    p_init     = Theta(:,1)';
    Theta(:,1) = [];
    b    = Theta(:,1)';
    %% Obtain the likelihood
    ret = 0;
    for i = 1:NC
        p = p_init;
        % prev_mean_est = mean(DATA(1,i,:))*p;
        Mij = get_Mean_TD_h(Theta,Conc(i),tv,Time);
        init_M = squeeze(DATA(1,i,:))*p;
        for j = 2:NT
            DATA_i = squeeze(DATA(j,i,:));  
            DATA_i = DATA_i(~isnan(DATA_i));
            NR     = length(DATA_i);
            % t       = Time(j) - Time(j-1);
            % Theta_t = Drug_T(Theta,tv,Time(j-1));
            % Aij      = Drug_A(Theta_t,Conc(i),cmd);
            % Mij      = get_Mean_alt(Aij,[0,t]); % s x s x NT tensor
            % Covariance 
            % Sigij   = zeros(s,s,s);
            % for k = 1:s
            %     temp = get_sig_j(Aij,b,t,k);
            %     Sigij(:,:,k) = temp;
            % end
            % prev_est = [];
            for r = 1:NR
                init = init_M(r,:);
                Meani = sum(init*Mij(:,:,j-1));
                % prev_est = [prev_est;init*Mij];
                % Vari   = get_Var_alt(init,Sigij);
                Var_ob_i =  c^2;
                Var_ob_i_inv = Var_ob_i^(-1);
                ret = ret + 1/2 * (DATA_i(r) - Meani)^2*Var_ob_i_inv + 1/2*log(2 * pi * Var_ob_i);
            end
            % prev_mean_est = mean(prev_est);
            % p = prev_mean_est./sum(prev_mean_est);
        end
    end