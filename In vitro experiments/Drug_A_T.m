function ret = Drug_A_T(theta,dose,t,tv)
    s = size(theta,1);
    if size(theta,2) ~= s+5
        warning('The dimension is not matched.')
        pause
    end
    nu  = theta(:,3:s+1);
    alpha_vec = theta(:,1);
    beta_vec  = theta(:,2);
    b_beta_vec = theta(:,s+2);
    E_beta_vec = theta(:,s+3);
    b_nu_vec   = theta(:,s+4);
    E_nu_vec   = theta(:,s+5);
    k  = tv(1);
    t0 = tv(2);
    rate = 1/(1+exp(tv(1)*(tv(2)-t)));
    b_beta_vec = 1 - rate + rate*b_beta_vec;
    b_nu_vec = 1 - rate + rate*b_nu_vec;
    for i = 1:s
        hill_beta_i = get_Hill2(dose,b_beta_vec(i),E_beta_vec(i));
        beta_vec(i) = beta_vec(i)-log(hill_beta_i);
    end
    lambda_vec = alpha_vec - beta_vec;
    ret   = diag(lambda_vec);
    for i = 1:s
        hill_nu_i = get_Hill2(dose,b_nu_vec(i),E_nu_vec(i));
        nu(i,:) = max(nu(i,:) + log(hill_nu_i),0);
        ret(i,1:i-1) = nu(i,1:i-1);
        ret(i,i+1:end) = nu(i,i:end);
    end
end