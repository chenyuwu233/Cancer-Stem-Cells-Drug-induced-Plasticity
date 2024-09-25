%%  Setting the Concentration and Time points

init   = 1000;
s      = 2;
Conc = 10^(6)*[0  31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6)];
NT   = 13;
Time = [0 : NT-1]*3;
NC   = length(Conc);
NR   = 20;
cmd  = 'CSC_DIS';

% p_range      = [0.25,0.75];
beta_s_range      = [1e-3,0.9];
beta_d_range      = [1e-3,0.5];
lam_s_range       = [0,0.1];
lam_d_range       = [0,0.05];
nu_diff_range = [0.05,0.4];
b_beta_range    = [0.8,0.9];
b_nu_range      = [1+1e-6,1.1];
E_range         = [0.0625,2.5];
c_range         = [0,10];


c      = rand*10;
% p      = rand*(p_range(2)-p_range(1)) + p_range(1);
% alpha1 = rand*(alpha_range(2)-alpha_range(1)) + alpha_range(1);
% alpha2 = rand*(alpha_range(2)-alpha_range(1)) + alpha_range(1);
beta1  = rand*(beta_s_range(2)-beta_s_range(1)) + beta_s_range(1);
beta2  = rand*(beta_d_range(2)-beta_d_range(1)) + beta_d_range(1);
alpha1 = beta1 + rand*lam_s_range(2);
alpha2 = beta2 + rand*lam_d_range(2);
nu12   = rand*(nu_diff_range(2)-nu_diff_range(1)) + nu_diff_range(1);
nu21   = 0;
b1_beta = 1;
b2_beta = rand*(b_beta_range(2)-b_beta_range(1)) + b_beta_range(1);
b2_nu   = 2-b2_beta;
E2_beta = rand*(E_range(2)-E_range(1)) + E_range(1);
E2_nu = E2_beta;

theta   = [alpha1,beta1,nu12,1,1,1,1,...
              alpha2,beta2,nu21,b2_beta,E2_beta,b2_nu,E2_nu]

A = [alpha1-beta1,nu12;0,alpha2-beta2];


get_stable_p(A)


%% 

idx = 129;
if idx <= 60
    name = strcat('Results\In silico Base\Boot_CI_CSC_DIS_alt_',num2str(idx),'.mat');
else
    name = strcat('Results\In silico Base\PE_CSC_DIS_',num2str(idx),'.mat');
end
load(name)

%% Plot the drug effects

gr_beta = [];
gr_nu = [];
gr_all = [];

% c     = theta(end);
% theta(end) = [];
Theta = coeff_reshape(theta,cmd,'Gen',s);
% Theta_like = coeff_reshape(theta,cmd,'Like',s);

% Theta(2,6) = 1;
Theta(2,4) = 1;
A    = Drug_A(Theta,0,cmd);
p    = get_stable_p(A);
init_p = init*p;
Z = [];

for c = Conc
    % Theta_c = Drug_theta(Theta,c,cmd);
    A_c = Drug_A(Theta,c,cmd);
    gr = max(eig(A_c));
    gr_all = [gr_all,gr];
    M_c = get_Mean(A_c,Time);
    Zc = 1000;
    for i = 1:NT-1
        Zc = [Zc,sum(init_p*M_c(:,:,i))];
    end
    Z = [Z;Zc];
end




% for gr = gr_all
%     Z_gr = exp(gr.*Time);
%     Z = [Z;Z_gr];
% end
%%

GR1_color = [254 129 125];
GR2_color = [129 184 223];
GR1_color = GR1_color./255;
GR2_color = GR2_color./255;
Color     = [GR1_color;GR2_color];



%%

for i = 1:11
    plot3(Conc(i)*ones(1,13),Time,Z(i,:),'LineWidth',3,'Color',Color(2,:));
    hold on
end

ax = gca;
ax.XScale = 'log';

plot3(Conc,Time(13)*ones(1,11),Z(:,13),'-o','LineWidth',3,'Color',Color(1,:),'MarkerSize',10)

ax.XTick = Conc;
ax.YTick = Time;
xlabel('Concentration levels')
ylabel('Time')
zlabel('Cell count')
title('Cell growth patten under the drug')
ax.FontWeight = 'bold';
