%% Return the Covariance matrix for the CLT_BP.
%  The input information:
%  X_init: 1 x NS vector records the initial cell number for each sub-type
%  Time: 1 x NT time vector
%  lam:  1 x NS vector records the growth rate for each sub-type
%  r_nu: 1 x NS vector records the birth rate plus death rate for each
%  sub-type
%  cmd: commend indicator for approximated-exact approach

function ret = get_Var(X_init, Time, lam, r_nu,cmd)
    NS  = length(X_init); % number of sub group
    NT  = length(Time);   % number of time point

%         ret = zeros(NT-1);
%         for i = 1:NT-1 % Note we start from 2 because X(0) is deterministic. It has 0 covariance.
%             for j = 1:NT-1
%                 var = 0;
%                 for n = 1:NS
%                    for m = 2:min(i+1,j+1)
%                        t     = Time(m) - Time(m-1);
%                        var_n = r_nu(n) / lam(n) * (exp(2 * lam(n) * t) - exp(lam(n) * t));
%                        var   = var + X_init(n) * exp( lam(n) * (Time(i+1) + Time(j+1) - 2*Time(m))) * exp(lam(n) * Time(m-1)) * var_n;
%                    end
%                 end
%                 ret(i,j) = var;
%             end
%         end



  
    if cmd == 0
        ret = zeros(NT-1);
        for n = 1:NS
%             if lam(n) < 0 
%                 n
%                 lam(n)
%                 pause
%                 continue
%             end
            A   = zeros(NT-1);
            Sig = zeros(NT-1); 
            for i = 1:NT-1    
                t        = Time(i+1) - Time(i);
                Sig(i,i) = r_nu(n) / lam(n) * (exp(2 * lam(n) * t) - exp(lam(n) * t));
                for j = 1:i
                    A(i,j) = exp(lam(n) * (Time(i+1) - Time(j+1)) + lam(n) * Time(j)/ 2);
                end
            end
            ret = ret + X_init(n) * A * Sig * A';
        end
    else
        ret = zeros(NT);
        for n = 1:NS
            A   = zeros(NT-1);
            Sig = zeros(NT-1); 
            for i = 1:NT-1    
                t        = Time(i+1) - Time(i);
                Sig(i,i) = r_nu(n) / lam(n) * (exp(2 * lam(n) * t) - exp(lam(n) * t));
                for j = 1:i
                    A(i,j) = exp(lam(n) * (Time(i+1) - Time(j+1)) + lam(n) * Time(j)/ 2);
                end
            end
            if n == 1
                ret(1:NT-1, 1:NT-1) = ret(1:NT-1, 1:NT-1) + X_init(n) * A * Sig * A';
            elseif n == 2
                temp =  X_init(n) * A * Sig * A';
                ret(1:NT-2, 1:NT-2) = ret(1:NT-2,1:NT-2) + temp(1:NT-2, 1:NT-2);
                ret(1:NT-2,end) = ret(1:NT-2,end) + temp(1:NT-2,end);
                ret(end,1:NT-2) = ret(end,1:NT-2) + temp(end,1:NT-2);
                ret(end,end)    = ret(end,end) + temp(end,end);
            end
        end
    end
%     if cmd == 0
%         ret = ret + sig2 * eye(NT-1);
%     else
%         ret = ret + sig2 * eye(NT );
%     end
%     A
%     Sig
end