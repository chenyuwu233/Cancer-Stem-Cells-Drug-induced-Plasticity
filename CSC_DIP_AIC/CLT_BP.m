% This is function using the functional central limit theorem to
% approximate the branching process path probability. We have input X_ij:
% path data, theta_j: [{p}_n, {r,nu}_n], num_sub:number of sub-group
% X_ij: NT x 1 vector includes the total number of cell at each time point
% theta_j: 3NS x 1 vector includes the parameter: [{p}_n, {r,nu}_n]
% Time: 1 x NT vector includes the time points
% num_sub: number of sub-types
% sig2: observation noise
% cmd: special commend

function ret = CLT_BP(X_ij,theta_j,Time,num_sub,sig2,cmd)
    p      = theta_j(1:num_sub);
    X_init = round(X_ij(1) * p);
    lam    = [];
    r_nu   = [];
    NT     = length(Time);
    for i = 1:num_sub
        r_nu = [r_nu,theta_j(num_sub + 2*i - 1) + theta_j(num_sub + 2*i)];
        lam  = [lam,theta_j(num_sub + 2*i - 1) - theta_j(num_sub + 2*i)];
    end
    %% Compute the mean and co-variance given the time
    if cmd == 0
        Mean = X_init*exp(lam'*Time);% 1 x NT vector that record expected total number of cell at Time.
    else
        Mean = X_init*exp(lam'*Time(1:end-1));
        Mean = [Mean, X_init(1) * exp(lam(1) * Time(end)), X_init(2) * exp(lam(2) * Time(end))];
    end
    Var  = get_Var(X_init, Time, lam, r_nu,cmd);
    Mean = Mean(2:end);
    X    = X_ij(2:end);
    if cmd == 0
        Var  = Var + sig2^2 * eye(NT-1);
    else
        Var = Var + sig2^2 *eye(NT);
    end
%     Var  = Var;
    if det(Var) <=0 
        ret = 0;
    else
        ret  = det(2 * pi * Var)^(-1/2) * exp(-1/2 * (X - Mean')' * Var^(-1) * (X - Mean'));
    end
    if isinf(ret)||~isreal(ret)
        ret = 0;
    end
end