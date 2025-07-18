%% Description
%  Function to obtain the GR50 from corresponding b,E and maximum dosage.
%  One important assumption here is that n = 1 here
%  Hill equation: H(d) = log(b + (1-b)/(1+(d/E)^n))
%  
%  Input:
%  - b: coefficient related to the maximum drug effect
%  - E: coefficient related to the half effect dosage
%  - d: maximum dosage level we applied in the experiment
%
%  Output:
%  - ret: correspondent GR50 value.

function ret = GR50_h(b,E,d,n)
    r_m = ( log(b + (1 - b)/(1 + (d/E)^n)))/2;
    ret = E * ((exp(r_m ) - 1)/(b - exp(r_m)))^(1/n);
end