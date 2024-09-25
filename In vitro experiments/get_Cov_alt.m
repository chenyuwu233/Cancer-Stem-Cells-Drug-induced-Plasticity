 %% Description:
%  This is a function to obtain the elements of Covaraince matrix for
%  switching model with n sub-types
%
%  Required function:
%  get_sig_j.m
%  get_m.m
%
%  Input:
%  init: 1 x n vector that contains the initial value for each sub-types
%  A: n x n rate matrix.
%  t1: the first time point's index, t1 <= t2
%  t2: the second time point's index, allow the unequal interval.
%  Time: Time vector.
%  Mi: s x s x NT tensor that records all possible mean matrix
%
%  Output:
%  Cov(S(t1), S(t2))

function ret = get_Cov_alt(init,A,Sig,t1,t2,Time)
    if t1 == 1
        if t1 == t2
            ret = 1;
        else
            ret = 0;
        end
    else
        X0  = init;
        n   = length(X0);
        ret = 0;
        if t1 == t2
            Mi_12 = eye(n);
        else
            Mi_12 = get_Mean(A,[0,Time(t2) - Time(t1)]);
        end
        for i = 1:n
            ret  = ret + X0(i)*ones(1,n)*Sig(:,:,i,t1)*Mi_12*ones(n,1);
        end
    end
end