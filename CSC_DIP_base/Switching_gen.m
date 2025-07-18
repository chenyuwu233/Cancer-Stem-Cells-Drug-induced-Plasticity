%% Description:
% This is a generating function for the switching model.
%
% Inputs (suppose we have s sub-types)
%  - init: 1 x s vector that includes initial cell numbers
%  - theta: 1 x s*(d+1)+1 that includes the initial proportion, birth rate, death rate, and
%    switching rate among every sub-types, the last element is drug effect.
%                   [{p,alhpa,beta,{nu}_{s-1},b,E,(n)}_s,c]
%  - Time: 1 x NT vector that includes all the time points we collect the data.
%  - Conc: 1 x NC vector that includes the concentration points
%  - NR: number of replicates
%  - s: number of sub-types
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
%      - 'Hill2_death_2': Assume the death rate (beta)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently.
%
%      - 'Multi_Selective': Assume the switching rate (nu)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently. We applied the drug independently.
%
%      - 'Multi_Selective_d': Assume the death rate (beta)_i 'Hilly'
%      depends on the drug concentration with parameters (b,E)_i with two
%      drugs effect independently. We applied the drug independently.
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
%      - 'CSC_DIS': One drug which can increase death rate and switching rate of the DC.
%        In this experiment we have the following parameters set:
%        θ=[α_1,β_1,ν_12, α_2=0,β_2,ν_21,b_beta_2,E_beta_2,b_nu_2,E_nu_2,c]
%
% Output:
% ret: NT x NC x NR tensor that records the total number of cells at every
% time points




function ret = Switching_gen(init, theta, Time, Conc, NR, NC, NT, s,cmd)
    c     = theta(end);
    theta(end) = [];
    Theta = coeff_reshape(theta,cmd,'Gen',s);
    
    A    = Drug_A(Theta,0,cmd);
    p    = get_stable_p(A);



    init = init*p;
    %% Generating the data
    switch cmd
        case {'Linear_switching','Single_Broad','Single_Broad_d','Hill2_switching_death','CSC_DIS'}
            ret = zeros(NT,NC,NR);
            fprintf('Start 1 drug data generation\n')
            for i = 1:NC
                for j = 1:NR
                    Theta_i = Drug_theta(Theta,Conc(i),cmd);
                    path_i  = Switching_path(init,Theta_i,Time);
                    ret(:,i,j) = max(0,sum(path_i,1)+[0,round(normrnd(0,c,1,NT-1))]);
                end
            end
        case {'Hill2_switching_Hill2_death_i','Hill2_death_Hill2_death_i',...
                'Multi_Broad','Multi_Selective','Multi_Selective_d','Multi_All_nu'}
            ret = zeros(NT,NC,NR);
            fprintf('Start 2 drug(independent) data generation\n')
            for i = 1:NC
                for j = 1:NR
                    Theta_i = Drug_theta(Theta,Conc(i,:),cmd);
                    path_i  = Switching_path(init,Theta_i,Time);
                    ret(:,i,j) = max(0,sum(path_i,1)+[0,round(normrnd(0,c,1,NT-1))]);
                end
            end
        case {'Hill2_switching_Hill2_death','Hill2_switching_2','Hill2_death_2',...
                'Hill2_death_Hill2_death','Hill2_switching_Hill2_switching'}
            ret = zeros(NT,NC(1),NC(2),NR);
            fprintf('Start 2 drug(cross) data generation\n')
            for i = 1:NC(1)
                for j = 1:NC(2)
                    for k = 1:NR
                        Theta_ij = Drug_theta(Theta,[Conc(1,i),Conc(2,j)],cmd);
                        path_ij  = Switching_path(init,Theta_ij,Time);
                        ret(:,i,j,k) = max(0,sum(path_ij,1)+[0,round(normrnd(0,c,1,NT-1))]);
                    end
                end
            end

    end

    
end