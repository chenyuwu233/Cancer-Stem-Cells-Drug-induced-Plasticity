%% Description:
%  This is a function to obtain the Covaraince matrix for swithing model
%
%  Required function:
%  get_Cov.m
%  get_sig_j.m
%  get_m.m
%  
%  Input:
%  init: 1 x s vector of initial number of cells of each sub-type
%  A: s x s rate matrix
%  Sig: s x s x s tensor that store all the covariance.
%  Time: Time vector
%  Mi: s x s x NT tensor that records all possible mean matrix
%
%  Output:
%  Covariance matrix of [S(t1),...,S(tn)]


function ret = get_Cov_mat(init, A, Sig, Time, Mi)
    NT  = length(Time);
    ret = zeros(NT,NT);
    for i = 1:NT
        for j = 1:NT
%             fprintf('Covariance: %d, %d: \n', i,j)
            ret(i,j) = get_Cov(init,A,Sig,i,j,Time,Mi); 
        end
    end
end