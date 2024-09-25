%% Description:
%  Obtain the target b value for given effect and rate
%  
%  Inputs (suppose we have s sub-populations):
%  - effect: objective effect (This may come from major drug effect or the inherent rate of subpopulation)
%  - E_r: target E value
%  - dose: maximum dosage level
%  - rate: multiplier for the effect
%  
%  Outputs:
%  - ret: the desire b value


function ret = get_target_b(effect,E_r,dose,rate)
    target_Effect = exp(effect/rate);
    ret = ((1+dose/E_r)*target_Effect-1)*E_r/dose;
end