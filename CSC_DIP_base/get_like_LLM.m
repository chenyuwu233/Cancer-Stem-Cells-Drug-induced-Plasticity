%% Description:
%  This is a function to obtain the Covaraince matrix for swithing model
%  One important assumption: the time interval is identical.
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
%  - Time: 1 x NT vector that records all the time points
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


function ret = get_like_LLM(DATA,x,Time,Conc,NR,NC,NT,s,cmd)
    theta = x;
    c = theta(end)';
    ti = Time(2) - Time(1);
%     if ti ~= Time(3) - Time(2)
%         ti
%         pause
%     end

    theta(end) = [];
    Theta = coeff_reshape(theta,cmd,'Like',s);
    
    A    = Drug_A(Theta,0,cmd);
    p    = get_stable_p(A);
    


    %% Obtain the likelihood
    ret = 0;
    switch cmd
        case {'Linear_switching','Single_Broad','Hill2_switching_death','Single_Broad_d','CSC_DIS'}
            for i = 1:NC
                Ai      = Drug_A(Theta,Conc(i),cmd);
                % Mean
                Mi      = get_Mean(Ai,Time); % s x s x NT tensor
                MM      = [];
                for r = 1:NR
                    init = DATA(1,i,r)*p;
                    Meani   = zeros(1,NT-1);
                    for j = 1:NT-1
                        Meani(j) = sum(init*Mi(:,:,j));
                    end
                    MM = [MM;Meani];
                end
                Var_ob_i = c^2*eye(NT-1);
                Var_ob_i_inv = Var_ob_i^(-1);
                DATA_i = squeeze(DATA(2:end,i,:));
                ret = ret + 1/2 * trace((DATA_i' - MM)*Var_ob_i_inv*(DATA_i' - MM)');
                ret = ret + NR/2*log(det(2 * pi * Var_ob_i));
            end
            if ret < 0
                ret = 10^NR;
            end
            ret = real(ret);
    end
end