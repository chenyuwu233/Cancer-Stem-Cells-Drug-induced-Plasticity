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
%  - DATA: NT X NC X NR tensor that records all the data
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
%        θ=[α_1,β_1,ν_12, α_2=0,β_2,ν_21,b_beta_2,E_beta_2,b_nu_2,E_nu_2,c]
%
%  Output:
%  The Likelihood of the DATA


function ret = get_like_BaF3(DATA,x,Time,Conc,NR,NC,NT,s,cmd)
    theta = x;
    c = theta(end)';
    ti = Time(2) - Time(1);


    theta(end) = [];
    Theta = reshape(theta,[],s)';
    p     = Theta(:,1)';
    Theta(:,1) = [];
    b    = Theta(:,1)';

    %% Obtain the likelihood
    ret = 0;
    for i = 1:NC
        
        Ai      = Drug_A_h(Theta,Conc(i),cmd);
        % Mean
        Mi      = get_Mean(Ai,Time); % s x s x NT-1 tensor
        % Covariance 
        Sig_i   = zeros(s,s,s,NT);
        for j = 1:s
            for t = 1:NT
                temp = get_sig_j(Ai,b,Time(t),j);
                Sig_i(:,:,j,t) = temp;
            end
        end

        for r = 1:NR
            DATA_ij = squeeze(DATA(:,i,r));
            Trans_mat = eye(NT);
            Trans_mat(isnan(DATA_ij),:) = [];
            Time_ij = Time*Trans_mat';
            Sig_ij = Sig_i(:,:,:,~isnan(DATA_ij));
            NT_ij = length(Time_ij);
            if NT_ij <=1
                continue
            end
            
            if Time_ij(1) ~= 0
                Time_ij = Time_ij - Time_ij(1);
            end
            DATA_ij(isnan(DATA_ij)) = [];
            if isempty(DATA_ij)
                continue
            end
            

            init = DATA_ij(1)*p;
            DATA_ij(1) = [];
            M_ij = zeros(NT_ij-1,1);
            for j = 1:NT_ij-1
                tj = Time_ij(j+1);
                M_ij(j) = sum(init*Mi(:,:,tj/ti));
            end
            Vari   = get_Cov_mat_alt(init,Ai,Sig_ij,Time_ij);
            Vari   = Vari(2:end,2:end);
            Var_ob_i = Vari + c^2*eye(NT_ij-1);
            V_ij    = Var_ob_i;


            ret     = ret + 1/2*(DATA_ij - M_ij)'*V_ij^(-1)*(DATA_ij - M_ij) + 1/2*log(det(2 * pi * V_ij));
        end
        
    end
end