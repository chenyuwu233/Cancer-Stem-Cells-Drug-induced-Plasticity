%% Description:
%  This is a function to obtain the likelihood for statistical model
%  described in equation (23) of the manuscript.
%
%  Required function:
%  get_m.m
%  get_Mean_TD.m
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
%
%      - 'CSC_DIS': One drug which can increase death rate and switching rate of the DC.
%        In this experiment we have the following parameters set:
%        θ=[α_1,β_1,ν_12, α_2=0,β_2,ν_21,b_beta_2,E_beta_2,b_nu_2,E_nu_2,tv,c]
%
%  Output:
%  The Likelihood of the DATA


function ret = get_like_TD_meanflow(DATA,x,Time,Conc,NR,NC,NT,s,cmd)
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
        % Time-Time(1);
        Mij = get_Mean_TD(Theta,Conc(i),tv,Time-Time(1));
        % keyboard
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