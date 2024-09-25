%% Description:
%  Obtain the stable distribution from the rate matrix
%  
%  Inputs (suppose we have s sub-populations):
%  - A: s x s rate matrix, where diagonal is net-growth rate and
%  off-diagonal is the switching rate.
%  
%  Outputs:
%  - ret: 1 x s vector that indicate the theoritical stable distribution
%  of each sub-population

function ret = get_stable_p(A)
    [~,D,V] = eig(A);
    [~,lam_indi] = max(diag(D));
    ret = V(:,lam_indi)./sum(V(:,lam_indi));
    ret = ret';
end