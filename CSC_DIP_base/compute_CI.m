


ii = 85;

load_name = strcat('Result_GS/PE_CSC_DIS_',num2str(ii),'.mat');

load(load_name)

CI = [];

%%  Compute the Hessian



func_a = @(a) func([a,opt_xx_pe(2:end)]);
tic
[hess_a, ~] = hessdiag(func_a,opt_xx_pe(1));
t = toc

CI = [CI,[opt_xx_pe(1) + 1.96/sqrt(20*hess_a);opt_xx_pe(1) - 1.96/sqrt(20*hess_a)]];


func_a = @(a) func([a,opt_xx_pe(2:end)]);
tic
[hess_a, ~] = hessdiag(func_a,opt_xx_pe(1));
t = toc

CI = [CI,[opt_xx_pe(1) + 1.96/sqrt(20*hess_a);opt_xx_pe(1) - 1.96/sqrt(20*hess_a)]];


%%

theta_est = [theta(1:3),theta(8:9),theta(11:15)];
opt_est   = [opt_xx_pe(1:3),opt_xx_pe(8:9),opt_xx_pe(11:15)];
Hess_est  = diag(opt_Hessian)';
Hess_est  = [Hess_est(1:3),Hess_est(8:9),Hess_est(11:15)];

CI_est    = [opt_est + 1.96./sqrt(20.*Hess_est);
             opt_est - 1.96./sqrt(20.*Hess_est)];



%%
theta_est = [theta(1:3),theta(8:9),theta(11:14)];

RE_hist = [];
Hess_hist = [];

for ii = 1:100
    MLE = mean(CI_hist(:,:,ii));
    RE_hist = [RE_hist;abs(MLE - theta_est)./theta_est];
    Hess = (1./((CI_hist(1,:,ii) - MLE)./1.96)).^2;
    Hess_hist = [Hess_hist;Hess];
    CI_hist(1,:,ii) = MLE + 1.96./sqrt(Hess/init); 
    CI_hist(2,:,ii) = MLE - 1.96./sqrt(Hess/init);
end



%%



hit_rate = zeros(1,9);

for ii = 1:100

    hit = theta_est>=CI_hist(2,:,ii) & theta_est <= CI_hist(1,:,ii);
    hit_rate = hit_rate + hit;
end



%% Boot prctile CI hit

hit_rate = zeros(1,9);
ist      = 0;
for ii = 31:130
    load_name = strcat('Result_CI/Boot_CI_CSC_DIS_alt_',num2str(ii),'.mat');
    load(load_name)
    theta_est = [theta(1:3),theta(8:9),theta(11:14)];
    CI  = prctile(B_parameter,[2.5,97.5]);
    CI  = [CI(:,1:3),CI(:,8:9),CI(:,11:14)];
    hit = theta_est>=CI(1,:) & theta_est <= CI(2,:);
    if sum(hit) <9
        ist = ist+1;
        % keyboard
    end
    hit_rate = hit_rate + hit;
end


