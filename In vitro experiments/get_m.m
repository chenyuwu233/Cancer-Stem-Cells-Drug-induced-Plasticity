%% Description:
% This is a function to obtian the j-th row of the mean matrix. 
function ret = get_m(A,t,j)
    M_t = expm(t.*A);
    ret = M_t(j,:);
end