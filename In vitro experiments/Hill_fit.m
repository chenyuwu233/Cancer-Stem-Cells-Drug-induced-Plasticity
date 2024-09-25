function ret = Hill_fit(x,Conc,DATA)
    if length(Conc) ~= length(DATA)
        pause
    end
    ret = 0;
    for i = 1:length(Conc)
        ret = ret + (log(get_Hill2(Conc(i),x(1),x(2))) - DATA(i))^2;
    end
end