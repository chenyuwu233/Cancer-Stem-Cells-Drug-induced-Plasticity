%% Description:
%  Drug A dependency function designed for obtaining the drug affected theta.
%  
%  Inputs (suppose we have s sub-populations):
%  - theta: s x d matrix that includes all parameters
%  - tv:    variable to specify the effect dependence on the time
%  - t:     current time
%
%  Outputs:
%  - ret:   s x s matrix that includes all the drug affected
%  parameters.




function ret = Drug_T_h(theta, tv, t)
    if length(tv) == 1
        s = size(theta,1);
        rate = exp(t - tv)/(1+exp(t-tv));  % 1PL
        ret = zeros(size(theta));
        for i = 1:s
            temp = theta(i,:);
            temp(4) = 1 - rate + rate*temp(4);
            temp(7) = 1 - rate + rate*temp(7);
            ret(i,:) = temp;
        end
    elseif length(tv) == 2
        s = size(theta,1);
        rate = 1/(1+exp(tv(1)*(tv(2)-t)));  % 2PL
        ret = zeros(size(theta));
        for i = 1:s
            temp = theta(i,:);
            temp(4) = 1 - rate + rate*temp(4);
            temp(7) = 1 - rate + rate*temp(7);
            ret(i,:) = temp;
        end
    end
end