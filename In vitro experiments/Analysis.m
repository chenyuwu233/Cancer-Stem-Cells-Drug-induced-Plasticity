%% Check the expected structure


est = opt_xx_pe;
est(end) = [];

Theta = reshape(est,[],s)';
p     = Theta(:,1)';
Theta(:,1) = [];
b    = Theta(:,1)';
init_c = mean(DATA(1,1,:))*p;
SC_p = [];
TC_hist = [];
for i = 1:NC
    Ai = Drug_A_h(Theta,Conc(i),cmd);
    pi = [];
    tc = [];
    for j = 1:NT
        Mj = get_Mean(Ai,[0,Time(j)]);
        Ec = init_c*Mj;
        pi = [pi,Ec(1)/sum(Ec)];
        tc = [tc,sum(Ec)];
    end
    TC_hist = [TC_hist;tc];
    SC_p = [SC_p;pi];
end

%%
for i = potential
    if i ~= 0
        est = params_hist_pe(i,:);
        est(end) = [];
        
        Theta = reshape(est,[],s)';
        p     = Theta(:,1)';
        Theta(:,1) = [];
        b    = Theta(:,1)';
        init_c = init*p;
        SC_p = [];
        TC_hist = [];
        for i = 1:NC
            Ai = Drug_A_h(Theta,Conc(i),cmd);
            pi = [];
            tc = [];
            for j = 2:NT
                Mj = get_Mean(Ai,[0,Time(j)]);
                Ec = init_c*Mj;
                pi = [pi,Ec(1)/sum(Ec)];
                tc = [tc,sum(Ec)];
            end
            TC_hist = [TC_hist;tc];
            SC_p = [SC_p;pi];
        end
        SC_p
        est
        pause
    end
end



%%  Record the 
