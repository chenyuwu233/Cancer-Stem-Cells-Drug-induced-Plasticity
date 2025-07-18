

function [Mean,Cov] = get_Moment_DIP(Theta,init,Time,tv,b,dose,NR,cmd)
    %%
    if Time(1) == 0
        Time = Time(2:end);
    end

    if tv(1) ==0
        A = Drug_A(Theta,dose,cmd);

        M = @(t) expm(t*A);
    else
        A_t = @(t) Drug_A_T(Theta,dose,t,tv);
        M = @(t) expm(integral(A_t,0,t,"ArrayValued",true));
    end
    %%
    if size(init,1) == 1

        Mean = zeros(1,length(Time));
        for i = 1:length(Time)
            Mean(i) = sum(init*M(Time(i)));
        end
        Mean = repmat(Mean,NR,1);
    else
        Mean = zeros(NR,length(Time));
        for i = 1:length(Time)
            Mean(:,i) = sum(init*M(Time(i)),2);
        end
    end
    %%
    Cov = zeros(length(Time),length(Time));
end