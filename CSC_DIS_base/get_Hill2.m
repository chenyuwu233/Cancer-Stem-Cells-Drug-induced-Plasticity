function ret = get_Hill2(dose,b,E)
    ret = b+(1-b)/(1+(dose/E));
end