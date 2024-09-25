%% Description:
%  This is a function that reorder the parameter under the stable p
%  assumption according to every sub-populations' drug effect.
%  
%  Input (suppose we have s sub-types):
%  - theta: 1 x d vector that record all the parameter about the
%  sub-population
%  - cmd1: model command
%  - cmd2: usage command
%  - s: number of sub-populations
%
%  Output:
%  - ret: desired theta

function ret = coeff_reorder(theta, cmd1,cmd2,s)
    Theta = coeff_reshape(theta,cmd1,cmd2,s);
    switch cmd1
        case {'Single_Broad','Hill2_switching_2'}
            Theta = sortrows(Theta,5);
    end
    ret   = [];
    for i = 1:s
        ret = [ret,Theta(i,:)];
    end
end