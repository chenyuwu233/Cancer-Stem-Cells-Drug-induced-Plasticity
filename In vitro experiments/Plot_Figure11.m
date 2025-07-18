%%

Conc = [0,0.125,0.25,0.5,1,2,4,8,16];
Time = [0,12,36];
NT   = length(Time);
NC   = length(Conc);
NR   = 2;

Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1F');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 1E');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 1E');

DATA_SC = [cell2mat(Cell_SC(2:10,2:3))];
DATA_TC = [cell2mat(Cell_TC(2:10,2:3))];
DMSO_1  = cell2mat(DMSO_TC(2,6:10));
DMSO_2  = cell2mat(DMSO_TC(2,11:15));


% 
DATA    = zeros(NT,NC,NR);
% DATA(1,:,1) = mean(DMSO_1)/2.4655; % [0-48] rate = 0.0188
% DATA(1,:,2) = mean(DMSO_2)/2.4655;
% DATA(1,:,1) = mean(DMSO_1)/1.9676; % [12-48]
% DATA(1,:,2) = mean(DMSO_2)/1.9676;
% DATA(1,:,1) = mean(DMSO_1)/1.5702; % [24-48]
% DATA(1,:,2) = mean(DMSO_2)/1.5702;

% load('rate24.mat')
% DATA(1,:,1) = mean(DMSO_1)./rate;
% DATA(1,:,2) = mean(DMSO_2)./rate;

DATA(3,:,1) = DATA_TC(:,1)*mean(DMSO_1)/100;
DATA(3,:,2) = DATA_TC(:,2)*mean(DMSO_2)/100;


load('Forward_D24_v.mat')
DATA(2,:,1) = DATA_rate(1,:);
DATA(2,:,2) = DATA_rate(2,:);

load('T12_2T_v.mat')
DATA(1,:,1) = T12_2T(1);
DATA(1,:,2) = T12_2T(2);

DATA = round(DATA);

Time = Time + 12;


%%

figure

hold on
lines = [];
for i = 1:9
    scatter(Time(1:2),mean(DATA(1:2,i,:),3),500,'blue','x','LineWidth',3)
    scatter(Time(3),mean(DATA(3,i,:),3),300,'black','o','LineWidth',3);
    line = plot(Time,mean(DATA(:,i,:),3),'LineWidth',4);
    lines = [lines,line];

end

ax = gca;

Conc_name = {'0 \mu M','0.125 \mu M','0.25 \mu M','0.5 \mu M','1 \mu M','2 \mu M','4 \mu M','8 \mu M','16 \mu M'};
xlabel('Time (hours)')
ylabel('Cell count')

ax.FontSize = 30;
ax.FontWeight = 'bold';
legend(lines,Conc_name,'Location','northwest');
