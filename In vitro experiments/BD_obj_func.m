%% Description
%  This is the function that compute the negative log likelihood of Mixture data
%  by using the CLT approximated.
%  DATA:     with the dimension [Time, Dosage, Replicant]
%  Time:     record all the time point for the data.
%  Dosage:   record all the dosage that we use for differnet replicant
%  theta:    [(p,r,nu,b,E,n)_i,sig] 
%  num_sub:  number of sub population
%  sig:      observation noise ~ N(0,sig)
%  We will return the log-likelihood function for the data set.

function ret = BD_obj_func(DATA, Conc, Time, NR, theta,num_sub)
    NC       = length(Conc);

    if num_sub > 1
        p_vec    = theta(1:6:6*num_sub-6);
        p_vec    = [p_vec',1-sum(p_vec)];
    else
        p_vec    = 1;
    end
    ret   = 0;
    idx = 0;
    for j = 1:NC
        d = Conc(j);
        theta_j = p_vec;
        for s = 1:num_sub % Computing the birth rate and death rate under drug d.
           if s ~=num_sub
               r_s  = theta(6*s-4);
               nu_s = theta(6*s-3);  
               b_s  = theta(6*s-2);
               E_s  = theta(6*s-1);
               n_s  = theta(6*s);
               nu_s = nu_s - log(b_s + (1 - b_s)/(1 + (d/ E_s)^n_s));
           else
               r_s  = theta(6*s-5);
               nu_s = theta(6*s-4);  
               b_s  = theta(6*s-3);
               E_s  = theta(6*s-2);
               n_s  = theta(6*s-1);
               nu_s = nu_s - log(b_s + (1 - b_s)/(1 + (d/ E_s)^n_s));
           end
           theta_j = [theta_j,[r_s,nu_s]];
        end
        sig2 = theta(6*num_sub);
        for i = 1:NR     
            X_ij = squeeze(DATA(:,j,i));
            indi = isnan(X_ij);
            X_ij(indi) = [];
            if isempty(X_ij) 
                continue
            end
            X_ij
            idx = idx + 1;
            Tij  = Time;
            Tij(indi) = [];
            p_ij = CLT_BP(X_ij,theta_j,Tij,num_sub,sig2,0);
            
            if p_ij ~=0
                ret = ret - log(p_ij);
            else
                ret = ret + 1e6;
            end
        end
    end
    if ret == 0
        ret = 1e6;
    end
    idx
%     ret_theta = Theta;
end