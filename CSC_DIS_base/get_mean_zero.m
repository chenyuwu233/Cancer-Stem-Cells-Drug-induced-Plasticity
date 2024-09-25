function ret = get_mean_zero(A,t)
    ret = [exp(A(1,1)*t), A(1,2)/(A(1,1)-A(2,2))*(exp(A(1,1)*t)-exp(A(2,2)*t));0,exp(A(2,2)*t)];
end