%% Description:
%  Drug dependency function designed for obtaining the drug affected theta.
%  
%  Inputs (suppose we have s sub-populations):
%  - theta: s x d matrix that includes all parameters
%  - dose:  dosage level

%  Outputs:
%  - ret:   s x (s+1) matrix that includes all the drug affected
%  parameters.


function ret = Drug_theta_hierarchy(theta,dose,cmd)
    s = size(theta,1);
    if size(theta,2) ~= s+5
        warning('The dimension is not matched.')
        pause
    end
    ret = theta(:,1:s+1);
    b_beta_vec = theta(:,s+2);
    E_beta_vec = theta(:,s+3);
    b_nu_vec   = theta(:,s+4);
    E_nu_vec   = theta(:,s+5);
    for i = 1:s
        hill_beta_i  = get_Hill2(dose,b_beta_vec(i),E_beta_vec(i));
        hill_nu_i    = get_Hill2(dose,b_nu_vec(i),E_nu_vec(i));
        ret(i,2)     = ret(i,2)-log(hill_beta_i);
        ret(i,3) = max(ret(i,3)+log(hill_nu_i),0);
    end

end