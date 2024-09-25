function Hill_est = Fit_hill_3(Conc, Drug_effect)


    
    func = @(x) Hill3_fit(x,Conc,Drug_effect);
    
    lb = [0.5,0,0.2];
    ub = [1,16,5];
    num_optim = 20;
    x_init = [];
    for i = 1:num_optim
        xi = rand(1,length(ub)).*(ub-lb)+lb;
        x_init = [x_init;xi];
    end
    
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','iter','algorithm','sqp');
    
    fval_hist_pe = [];
    params_hist_pe = [];
    
    for j = 1:num_optim
    
        [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),[],[],[],[],lb,ub,[],options1); 
        fval_hist_pe = [fval_hist_pe,ff];
        params_hist_pe = [params_hist_pe;xx];
    
    end
    
    [of,oi] = min(fval_hist_pe);
    of
    Hill_est = params_hist_pe(oi,:);

end