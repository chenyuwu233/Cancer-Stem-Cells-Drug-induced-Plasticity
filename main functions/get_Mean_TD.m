%% Description:
%  This is a function to obtain all possible mean matrix used.
%
%  Assumption: the Time intervals are identical.
%
%  Input:
%  Theta: 2 x 7 matrix, each row corresponding to parameters of
%  one subpopulation.
%  dose: dosage applied.
%  tv: logistic time-delayed variables k, t0.
%  Time: 1 x NT Time vector.
%
%  Output:
%  s x s x NT-1 mean matrix for NT-1 intervals.
%
% ************ This is only for 2 subpopulation situation **************




function ret = get_Mean_TD(Theta,dose,tv,Time)
    NT  = length(Time);
    n   = size(Theta,1);
    ret = zeros(n,n,NT-1);

    
    for i = 2:NT
        A_t = @(t) Drug_A_T(Theta,dose,t,tv);
        A  = integral(A_t,0,Time(i),"ArrayValued",true);
        Mi = expm(A);
        ret(:,:,i-1) = Mi;
    end
end