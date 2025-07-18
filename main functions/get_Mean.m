%% Description:
%  This is a function to obtain all possible mean matrix used.
%
%  Assumption: the Time intervals are identical.
%
%  Input:
%  A: s x s rate matrix
%  Time: Time vector
%
%  Output:
%  s x s x NT-1 mean matrix for NT-1 intervals.




function ret = get_Mean(A,Time)
    NT  = length(Time);
    ret = zeros(size(A,1),size(A,2),NT-1);
    
    for i = 2:NT
         Mi = expm((Time(i) - Time(1)).*A);
         ret(:,:,i-1) = Mi;
    end
end