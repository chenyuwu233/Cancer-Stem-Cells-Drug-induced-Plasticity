%% Function for random event searching grid
% This is function that find the specific grid for random event simulation
%
% Inputs 
% rd_event: random number belongs to range that represents the uniform random
% event.
% range: Target range.
% num_vec: cell number of each grid, which can also be consider as weight.
% rate: rate for each individual cell or each unit of weight.
%
% Output:
% ret: index that the random event belongs to.

function ret = search_grid(rd_event,range,num_vec)
    if rd_event < range(1) || rd_event > range(2)
        rd_event
        range
        pause
    end
    rate = (range(2)-range(1))/sum(num_vec);
    cum_vec = range(1)+cumsum(num_vec.*rate);
    if rd_event < cum_vec(1)
        ret = 1;
    else
        idx = (cum_vec< rd_event);
        ret = find(idx,1,'last')+1;
        if ret > length(num_vec)
            cum_vec
            range
            pause
        end
    end
end