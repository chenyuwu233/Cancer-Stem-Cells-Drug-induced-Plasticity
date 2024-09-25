%% Energy distance
%  Compute the energy distance between two distribution F(x) and G(x)
%  Input:
%  D1: n x p matrix that records the data sample from the first distribution
%  F(x)
%  D2: m x p matrix that records the data sample from the second
%  distribution G(X)
% 
%  Output:
%  ret: The energy distance between two distribution
%
%

function ret = Energy_dist(D1,D2)
    if size(D1,2) ~= size(D2,2)
        error('The dimension of two distribution are not the same')
    end
    n = size(D1,1);
    m = size(D2,1);
    ret = 0;
    %% A
    for i = 1:n
        for j = 1:m
            ret = ret + 2/(n*m)*norm(D1(i,:) - D2(j,:));
        end
    end
    %% B
    for i = 1:n
        for j = 1:n
            ret = ret - 1/(n^2)*norm(D1(i,:) - D1(j,:));
        end
    end
    %% C
    for i = 1:m
        for j = 1:m
            ret = ret - 1/(m^2)*norm(D2(i,:) - D2(j,:));
        end
    end
    ret = sqrt(ret);
end