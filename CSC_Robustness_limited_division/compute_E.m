%% Description:
% This is a function to compute the expected number of descendents from a
% single first generation differentiated cell.
%
% Inputs (suppose we have 2 sub-types)
% b: birth rate of the differentiated cell.
% d: death rate of the differentiated cell.
% nu: de-differentiated rate.
% T: target time points
% max_cap: maximum descendent of the first generation differentiated cell.
% cmd: Whether this is for switching or asymmetric birth.
%
% Output:
% ret: expected number of descendents froma single first generation
% differentiated cell.

function ret = compute_E(b,d,nu,T,max_cap,cmd)
    switch cmd
        case 'switching'
            C   = (2*b/(nu+b));
            ret = (C^max_cap)*exp(-d*T); % C1 exp(-beta t)
            if max_cap > 0
                gamma = b + d + nu;
                ret = ret + (1-C^max_cap)*exp(-gamma*T); % C2 exp(-gamma t)
                for i = 1:max_cap-1
                    ret = ret + (2*b*T)^i*(1/factorial(i))*(1-C^(max_cap-i))*exp(-gamma*T);
                end
            end
        case 'asymmetric'
            C   = 2;
            ret = (C^max_cap)*exp(-d*T); % C1 exp(-beta t)
            if max_cap > 0
                gamma = b + d;
                ret = ret + (1-C^max_cap)*exp(-gamma*T); % C2 exp(-gamma t)
                for i = 1:max_cap-1
                    ret = ret + (2*b*T)^i*(1/factorial(i))*(1-C^(max_cap-i))*exp(-gamma*T);
                end
            end
    end
end