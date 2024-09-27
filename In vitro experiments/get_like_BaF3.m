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