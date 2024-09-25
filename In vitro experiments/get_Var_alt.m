%% Return the Covariance matrix for the CLT_BP.
%  This is a function to obtain the variance for swithing end-points model
%
%  Required function:
%  
%  
%  Input:
%  init: 1 x s vector of initial number of cells of each sub-type
%  Sig: s x s x s tensor that store all the covariance.
%  vec: 1 x s vector which indicate the specific subpopulation variance we
%  are interested in obtaining.
%
%  Output:
%  Variance of S(t_i)

function ret = get_Var(init,Sig,vec)
    s = size(Sig,3);
    ret = 0;
    if nargin == 3
        for i = 1:s
            ret = ret + init(i)*vec'*Sig(:,:,i)*vec;
        end
    else
        for i = 1:s
            ret = ret + init(i)*ones(1,s)*Sig(:,:,i)*ones(s,1);
        end
    end
end