%% Function for random event searching grid
% This is function that output a random index
%
% Inputs 
% vec: vector that we would like to find a random index that larger than 0
%
% Output:
% ret: index that the random event belongs to.

function ret = get_rand_index(vec)
    idx = randi(sum(vec>0));
    vec_idx = find((vec>0),idx);
    ret = vec_idx(end);
end