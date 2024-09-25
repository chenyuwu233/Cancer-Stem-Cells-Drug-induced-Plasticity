%% Description:
% This is a function that return the matrix D when computing the covariance
% matrix of multi-type branching process generated from single type j cell at the time t.
%
% input:
% - A: rate matrix
% - t: Time variable
% - b: Birth rate vector for every sub-types
% - j: Target sub-type

function ret = get_sig_j(A,b,t,j)
    Ab  = A - diag(diag(A)) + diag(b);
    M   = @(x) expm(x .* A);
    m_j = get_m(A,t,j);
    obj = @(tau) M(t - tau)' * (Ab.*get_m(A,tau,j)' + (Ab.*get_m(A,tau,j)')') * M(t - tau);
    D   =  integral(obj,0,t,"ArrayValued",true);
    ret = D + diag(m_j) - (m_j)'*m_j;
end