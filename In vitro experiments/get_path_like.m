%% Description:
%  This is a function to obtain the Covaraince matrix for swithing model
%
%  Required function:
%  get_Cov_mat.m
%  get_Cov.m
%  get_sig_j.m
%  get_m.m
%  
%  Input:
%  path: 1 x NT vector that records the path 
%  A: s x s rate matrix
%  Sig: s x s x s tensor that store all the covariance.
%  Time: Time vector
%  init: 1 x s vector that records the initial cell counts
%  Mi: s x s x NT vector that represent the mean alone the time.
%
%  Output:
%  The probability of obtaining the given path


function ret = get_path_like(path,A,Mi,Sig,Time,init,ob_std)
    X    = path;
    NT   = length(Time);
    Var  = get_Cov_mat(init,A,Sig,Time,Mi);
    Var  = Var(2:end,2:end);
    Var_ob = Var + ob_std^2*eye(NT-1);
    Mean = zeros(1,NT-1);
    for i = 1:NT-1
        Mean(i) = sum(init*Mi(:,:,i));
    end  
    X    = X(2:end);
    ret  = -1/2*log(det(2 * pi * Var_ob)) - 1/2 * (X - Mean) * Var_ob^(-1) * (X - Mean)';
end