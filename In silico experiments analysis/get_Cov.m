 %% Description:
%  This is a function to obtain the elements of Covaraince matrix for
%  switching model with n sub-types
%
%  Required function:
%  get_sig_j.m
%  get_m.m
%
%  Input:
%  init: 1 x n vector that contains the initial value for each sub-types
%  A: n x n rate matrix.
%  t1: the first time point's index
%  t2: the second time point's index
%  Time: Time vector.
%  Mi: s x s x NT tensor that records all possible mean matrix
%
%  Output:
%  Cov(S(t1), S(t2))

function ret = get_Cov(init,A,Sig,t1,t2,Time,Mi)
    t   = min(t1,t2);
    ti  = Time(2) - Time(1);
    if t == 1
        if t1 == t2
            ret = 1;
        else
            ret = 0;
        end
    else
        X0  = init;
        n   = length(X0);
        ret = 0;
    %     for i = 1:n
    %         temp = 0;
    %         for j = 2:t 
    %             for k = 1:n
    %                 tm = Time(j-1);    % t_(m-1)
    %                 t1m = Time(t) - Time(j); % t_L - t_m
    %                 tmm = Time(j) - Time(j-1);
    %                 M1 = expm(tm .*A);
    %                 M2 = expm(t1m .*A);
    %                 S1 = get_sig_j(A,b,tmm,k);
    %                 temp = temp + M1(i,k) * (ones(1,n) * M2' * S1 * M2 * ones(n,1));
    %             end
    %         end
    %         temp
    %         pause
    %         ret = ret + X0(i) * temp;
    %     end
        for i = 1:n
            temp = 0;
            for m = 2:t 
                for k = 1:n
                    tm = Time(m-1);    % t_(m-1)
                    t1m = Time(t1) - Time(m); % t_1 - t_m
                    t2m = Time(t2) - Time(m); % t_2 - t_m
%                     tmm = Time(m) - Time(m-1); % t_m - t_(m-1)
%                     M1 = expm(tm .*A);
%                     M2 = expm(t1m .*A);
%                     M3 = expm(t2m .*A);
                    
                    if tm == 0
                        M1 = eye(size(A));
                    else
                        M1 = Mi(:,:,tm/ti);
                    end
                    
                    if t1m == 0
                        M2 = eye(size(A));
                    else
                        M2 = Mi(:,:,t1m/ti);
                    end
                    
                    if t2m ==0
                        M3 = eye(size(A));
                    else
                        M3 = Mi(:,:,t2m/ti);
                    end
%                     S1 = get_sig_j(A,Sig,tmm,k);
%                     S1
%                     pause
                    S1 = Sig(:,:,k);
%                     eig(S1)
                    temp = temp + M1(i,k) * trace(M2 * ones(n,1) * ones(1,n) * M3' * S1);
                end
            end
%             fprintf('Type %d temp: \n', i)
%             X0(i) * temp
%             pause
            ret = ret +  X0(i)*temp;
        end
    end
end