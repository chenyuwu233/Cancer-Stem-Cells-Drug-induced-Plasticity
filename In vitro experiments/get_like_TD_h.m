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
%
%  Output:
%  The Likelihood of the DATA


function ret = get_like_TD_h(DATA,x,Time,Conc,NR,NC,NT,s,cmd)
    theta = x;
    c = theta(end)';
    

    theta(end) = [];
    if length(theta) == 22
        tv = theta(end-1:end); % Variable to control the effect delay
        theta(end-1:end) = [];
    else
        tv = theta(end); % Variable to control the effect delay
        theta(end) = [];
    end


    Theta = reshape(theta,[],s)';
    p     = Theta(:,1)';
    Theta(:,1) = [];
    b    = Theta(:,1)';
    
    % init = DATA(1,1,1)*p; % Assume the unify initial cell count
    %% Obtain the likelihood
    ret = 0;
    for i = 1:NC
        for j = 2:NT
            DATA_i = squeeze(DATA(j,i,:));  
            DATA_i = DATA_i(~isnan(DATA_i));
            NR     = length(DATA_i);
            t       = Time(j);
            Theta_t = Drug_T_h(Theta,tv,t);
            
            Aij      = Drug_A_h(Theta_t,Conc(i),cmd);
            % Mean
            Mij      = get_Mean_alt(Aij,[0,t]); % s x s x NT tensor
            % Covariance 
            Sigij   = zeros(s,s,s);
            for k = 1:s
                temp = get_sig_j(Aij,b,t,k);
                Sigij(:,:,k) = temp;
            end

            for r = 1:NR
                init = DATA(1,i,r)*p;
                Meani = sum(init*Mij);
                Vari   = get_Var_alt(init,Sigij);
                Var_ob_i = Vari + c^2;
                Var_ob_i_inv = Var_ob_i^(-1);
                ret = ret + 1/2 * (DATA_i(r) - Meani)^2*Var_ob_i_inv + 1/2*log(2 * pi * Var_ob_i);
            end
        end
    end