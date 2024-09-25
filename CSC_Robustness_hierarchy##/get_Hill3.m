function ret = get_Hill3(dose,b,E,n)
    ret = b+(1-b)/(1+(dose/E)^n);
end