%% Description:
% This is a generating function for the switching model. 
% Special assumption: Hierarchy and initial p is all stem cells
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
%
%
% Output:
% ret: NT x NC x NR tensor that records the total number of cells at every
% time points




function ret = Switching_gen_hierarchy(init, theta, Time, Conc, NR, NC, NT, s,cmd,varargin)
    input_num = nargin;
    if input_num == 10
        max_cap = varargin{1};
    else
        max_cap = 1;
    end
    s     = s + max_cap;
    p = zeros(1,s);
    p(1) = 1;

    theta = hierarchy_translate(theta,max_cap);
    c     = theta(end);
    theta(end) = [];
    
    Theta = reshape(theta,[],s)';


    init = round(init*p);
    %% Generating the data
    ret = zeros(NT,NC,NR);
    fprintf('Start 1 drug data generation\n')
    for i = 1:NC
        ci = Conc(i);
        for j = 1:NR
            Theta_i = Drug_theta(Theta,ci,cmd);
%                     Theta_i
%                     pause
%                     tic
            
            path_i  = Switching_path(init,Theta_i,Time);
%                     tj = toc
            ret(:,i,j) = max(0,sum(path_i,1)+[0,round(normrnd(0,c,1,NT-1))]);
        end
    end
end