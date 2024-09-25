function ret = get_Hill4(dose,d,b,E,n)
    ret = b+(d-b)/(1+(dose/E)^n);
end