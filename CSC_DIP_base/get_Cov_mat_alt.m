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


function ret = get_Cov_mat_alt(init,A,Sig,Time)
    NT  = length(Time);
    ret = zeros(NT,NT);
    for i = 1:NT
        for j = i:NT
            ret(i,j) = get_Cov_alt(init,A,Sig,i,j,Time); 
            if j ~= i
                ret(j,i) = ret(i,j);
            end
        end
    end
end