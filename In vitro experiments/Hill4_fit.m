function ret = Hill4_fit(x,Conc,DATA)
    if length(Conc) ~= length(DATA)
        pause
    end
    ret = 0;
    for i = 1:length(Conc)
        ret = ret + (get_Hill4(Conc(i),x(1),x(2),x(3),x(4)) - DATA(i))^2;
    end
end