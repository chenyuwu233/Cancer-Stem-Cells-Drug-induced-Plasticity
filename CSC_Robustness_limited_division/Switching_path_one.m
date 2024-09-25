%% Description:
% This is a path generating function for the switching model.
%
% Inputs (suppose we have n sub-types)
% init: 1 x n vector that includes initial cell numbers
% theta: n x n+1 matrix that includes the birth rate, death rate, and
% switching rate among every sub-types
% Time: 1 x t vector that includes all the time points we collect the data.
%
% Output:
% ret: n x t vector that records the number of each type cells at every
% time points



function ret = Switching_path_one(init,theta,Time)
    n = length(init);
    X = init'; % n x 1 vector that represents the current cell number.
    if size(theta,1) ~= n
        size(theta,1)
        pause
    end
    a = []; % 1 x n vector of the rate
    for i = 1:n
        ai = sum(theta(i,:));
        a = [a,ai];
    end
    t  = Time(1);
    ret = zeros(n,length(Time));
    ret(:,1) = init;
    nt = 2;
    id = 0;
    
    while t < Time(nt)&&sum(X) >0
        rate   = a * X;
        rate_v = a .* X';% 1 x n vector records the rate for each sub-types.
        t    = t - (1/rate) * log(rand);
        id = id + 1;
        if t > Time(nt)
            ret(:,nt) = X;
            nt      = nt+1;
            if nt > length(Time)
                break
            end
        end
        if sum(X) > 1e8
            ret(:,nt:end) = 1e8;
            break
        end
        event  = rand;
        indi   = 1:n;
        for i = 1:length(rate_v)
            block = sum(rate_v(1:i))/rate; 
            if event < block
                indi(i) = [];
                event_trans = rand;
                for j = 1:n+1
                    block_i = sum(theta(i,1:j))/sum(theta(i,:));
                    if event_trans < block_i
%                         t
%                         pause
                        if j == 1 % Birth
                            X(i) = X(i) + 1;
                        elseif j == 2 % Death
                            X(i) = X(i) - 1;
                        else % Birth of other sub-type cell
                            indi_j = indi(j - 2);
%                             X(i)   = X(i) - 1;
                            X(indi_j) = X(indi_j) + 1;
                        end
                        break
                    end
                end
                break
            end
        end
    end
end