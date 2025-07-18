%% Description:
%  Drug dependency function designed for obtaining the drug affected theta.
%  
%  Inputs (suppose we have s sub-populations):
%  - theta: s x d matrix that includes all parameters
%  - dose:  dosage level
%  - cmd:   string command about what model we used for drug dependency:
%      - 'Linear_switching': Assume the switching rate (nu)_i 'Linearly'
%      depends on the drug concentration with parameters (b)_i
%
%      - 'Single_Broad': Assume the switching rate (nu)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i
%
%      - 'Single_Broad_d': Assume the death rate (beta)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i
%
%      - 'Hill2_swtiching_2': Assume the switching rate (nu)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently.
%
%      - 'Multi_Selective':  Assume the switching rate (nu)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently.
%      Note that We applied the drug independently.
%
%      - 'Hill2_death_2': Assume the death rate (beta)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently.
%
%      - 'Multi_Selective_d': Assume the death rate (beta)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently.
%      Note that We applied the drug independently.
%
%      - 'Hill2_switching_death': Assume the switching rate (nu)_i and
%      death rate (beta)_i 'Hilly' depends on the drug concentration with
%      parameters (b_nu,E_nu,b_beta,E_beta)_i. 
%      Note that this model only has one type of drug.
%
%      - 'Hill2_switching_Hill2_death': Assume the switching rate (nu)_i
%      and death rate (beta)_i 'Hilly' depends on the drug concentration
%      with parameters (b_nu,E_nu,b_beta,E_beta)_i. 
%      Note that this model allows two type of drug, one affect the
%      switching rate (nu)_i and one affect the death rate (beta)_i.
%
%      - 'Hill2_switching_Hill2_switching': Assume the switching rate (nu)_i
%      'Hilly' depends on the drug concentration
%      with parameters (b_nu1,E_nu1,b_nu2,E_nu2)_i. 
%      Note that this model allows two type of drug, one affect the
%      switching rate with (b_nu1)_i and one affect the switching rate with
%      (b_nu2)_i
%
%      - 'Hill2_death_Hill2_death': Assume the death rate (beta)_i 
%      'Hilly' depends on the drug concentration
%      with parameters (b_beta1,E_beta1,b_beta2,E_beta2)_i. 
%      Note that this model allows two type of drug, one affect the
%      death rate with (b_beta1)_i and one affect the death rate with
%      (b_beta2)_i
%
%      - 'Hill2_switching_Hill2_death_i': Assume the switching rate (nu)_i
%      and death rate (beta)_i 'Hilly' depends on the drug concentration
%      with parameters (b_nu,E_nu,b_beta,E_beta)_i. 
%      Note that this model allows two type of drug, one affect the
%      switching rate (nu)_i and one affect the death rate (beta)_i. We
%      applied the drug independently.
%      
%      - 'Multi_Broad': Assume the switching rate (nu)_i
%      'Hilly' depends on the drug concentration
%      with parameters (b_nu1,E_nu1,b_nu2,E_nu2)_i. 
%      Note that this model allows two type of drug, one affect the
%      switching rate with (b_nu1)_i and one affect the switching rate with
%      (b_nu2)_i. We applied the drug independently.
%
%      - 'Hill2_death_Hill2_death_i': Assume the death rate (beta)_i 
%      'Hilly' depends on the drug concentration
%      with parameters (b_beta1,E_beta1,b_beta2,E_beta2)_i. 
%      Note that this model allows two type of drug, one affect the
%      death rate with (b_beta1)_i and one affect the death rate with
%      (b_beta2)_i. We applied the drug independently.      
%
%      - 'Multi_All_nu': Assume the death rate and asymmetrical birth rate
%      'Hilly' depends on the drug concentration but only significantly
%      affect nu of one subpopulation from one drug.
%
%      - 'SSDCC_d': One drug which can only eliminate the DCC.
%		 In this experiment we have the following parameters set:
%	     θ=[α_1,β_1,ν_12, α_2=0,β_2,ν_21,b_2,E_2,c]
%      
%      - 'CSC_DIS': One drug which can increase death rate and switching rate of the DC.
%        In this experiment we have the following parameters set:
%        θ=[α_1,β_1,ν_12, α_2=0,β_2,ν_21,b_beta_2,E_beta_2,b_nu_2,E_nu_2,c]

%
%  Outputs:
%  - ret:   s x (s+1) matrix that includes all the drug affected
%  parameters.


function ret = Drug_theta(theta,dose,cmd)
    s = size(theta,1);
    switch cmd 
        case {'Linear_switching','Linear_switching_EP'}
            if size(theta,2) ~= s+2
                warning('The dimension is not matched.')
                pause
            end
            ret = theta(:,1:s+1);
            b_vec = theta(:,s+2);
            for i = 1:s
                ret(i,3:end) = ret(i,3:end) + dose*b_vec(i);
            end
        case {'Single_Broad','Single_Broad_EP'}
            if size(theta,2) ~= s+3
                warning('The dimension is not matched.')
                pause
            end
            ret = theta(:,1:s+1);
            b_vec = theta(:,s+2);
            E_vec = theta(:,s+3);
            for i = 1:s
                hill_i = get_Hill2(dose,b_vec(i),E_vec(i));
                ret(i,3:end) = max(ret(i,3:end)+log(hill_i),0);
            end
        case {'Single_Broad_d','Hill2_death_EP'}
            if size(theta,2) ~= s+3
                warning('The dimension is not matched.')
                pause
            end
            ret = theta(:,1:s+1);
            b_vec = theta(:,s+2);
            E_vec = theta(:,s+3);
            for i = 1:s
                hill_i = get_Hill2(dose,b_vec(i),E_vec(i));
                ret(i,2) = ret(i,2)-log(hill_i);
            end
        case {'Hill2_switching_2','Multi_Selective','Hill2_switching_2_EP','Multi_Selective_EP'}
            if size(theta,2) ~= s+3
                warning('The dimension is not matched.')
                pause
            end
            if length(dose) ~= s
                warning('The drug dosage dimension is not enough')
                pause
            end
            ret = theta(:,1:s+1);
            b_vec = theta(:,s+2);
            E_vec = theta(:,s+3);
            for i = 1:s
                hill_i = get_Hill2(dose(i),b_vec(i),E_vec(i));
                ret(i,3:end) = max(ret(i,3:end)+log(hill_i),0);
            end
        case {'Hill2_death_2','Multi_Selective_d','Hill2_death_2_EP','Multi_Selective_d_EP'}
            if size(theta,2) ~= s+3
                warning('The dimension is not matched.')
                pause
            end
            if length(dose) ~= s
                warning('The drug dosage dimension is not enough')
                pause
            end
            ret = theta(:,1:s+1);
            b_vec = theta(:,s+2);
            E_vec = theta(:,s+3);
            for i = 1:s
                hill_i = get_Hill2(dose(i),b_vec(i),E_vec(i));
                ret(i,2) = ret(i,2)-log(hill_i);
            end
        case {'Hill2_switching_death','Hill2_switching_death_EP','CSC_DIS'}
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
                ret(i,3:end) = max(ret(i,3:end)+log(hill_nu_i),0);
            end
        case {'Hill2_switching_Hill2_death','Hill2_switching_Hill2_death_i','Hill2_switching_Hill2_death_EP','Hill2_switching_Hill2_death_i_EP'}
            if size(theta,2) ~= s+5
                warning('The dimension is not matched.')
                pause
            end
            if length(dose) ~= 2
                warning('The dimension of dosage level is not matched.')
                pause
            end
            b_beta_vec = theta(:,s+2);
            E_beta_vec = theta(:,s+3);
            b_nu_vec   = theta(:,s+4);
            E_nu_vec   = theta(:,s+5);
            ret = theta(:,1:s+1);
            for i = 1:s
                hill_beta_i  = get_Hill2(dose(1),b_beta_vec(i),E_beta_vec(i));
                hill_nu_i    = get_Hill2(dose(2),b_nu_vec(i),E_nu_vec(i));
                ret(i,2)     = ret(i,2)-log(hill_beta_i);
                ret(i,3:end) = max(ret(i,3:end)+log(hill_nu_i),0);
            end
        case {'Hill2_switching_Hill2_switching','Multi_Broad','Hill2_switching_Hill2_switching_EP','Multi_Broad_EP'}
            if size(theta,2) ~= s+5
                warning('The dimension is not matched.')
                pause
            end
            if length(dose) ~= 2
                warning('The dimension of dosage level is not matched.')
                pause
            end
            b_nu1_vec = theta(:,s+2);
            E_nu1_vec = theta(:,s+3);
            b_nu2_vec   = theta(:,s+4);
            E_nu2_vec   = theta(:,s+5);
            ret = theta(:,1:s+1);
            for i = 1:s
                hill_nu1_i  = get_Hill2(dose(1),b_nu1_vec(i),E_nu1_vec(i));
                hill_nu2_i  = get_Hill2(dose(2),b_nu2_vec(i),E_nu2_vec(i));
                ret(i,3:end) = max(ret(i,3:end)+log(hill_nu1_i),0);
                ret(i,3:end) = max(ret(i,3:end)+log(hill_nu2_i),0);
            end
        case {'Hill2_death_Hill2_death','Hill2_death_Hill2_death_i','Hill2_death_Hill2_death_EP','Hill2_death_Hill2_death_i_EP'}
            if size(theta,2) ~= s+5
                warning('The dimension is not matched.')
                pause
            end
            if length(dose) ~= 2
                warning('The dimension of dosage level is not matched.')
                pause
            end
            b_beta1_vec = theta(:,s+2);
            E_beta1_vec = theta(:,s+3);
            b_beta2_vec   = theta(:,s+4);
            E_beta2_vec   = theta(:,s+5);
            ret = theta(:,1:s+1);
            for i = 1:s
                hill_beta1_i  = get_Hill2(dose(1),b_beta1_vec(i),E_beta1_vec(i));
                hill_beta2_i  = get_Hill2(dose(2),b_beta2_vec(i),E_beta2_vec(i));
                ret(i,2)     = ret(i,2)-log(hill_beta1_i);
                ret(i,2)     = ret(i,2)-log(hill_beta2_i);
            end
        case {'Multi_All_nu'}
            if size(theta,2) ~= s+3
                warning('The dimension is not matched.')
                pause
            end
            if length(dose) ~= s
                warning('The drug dosage dimension is not enough')
                pause
            end
            ret = theta(:,1:s+1);
            b_vec = theta(:,s+2);
            E_vec = theta(:,s+3);
            for i = 1:s
                hill_i = get_Hill2(dose(i),b_vec(i),E_vec(i));
                ret(i,3:end) = max(ret(i,3:end)+log(hill_i),0);
            end
    end
end