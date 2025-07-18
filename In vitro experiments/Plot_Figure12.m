COLO858_data = readcell('COLO858 data.xlsx','Sheet','COLO858-Vem live cell count'); % No provided ask if necessary.


%% Extract data to NT x NC x NR
Time = cell2mat(COLO858_data(3:end,1));
Conc = [0,0.032,0.1,0.32,1,3.2];
DATA_COLO858 = zeros(length(Time),length(Conc),4);
for i = 1:length(Conc)
    for j = 1:length(Time)
        DATA_COLO858(j,i,:) = cell2mat(COLO858_data(j+2,4*i-2:4*i+1));
    end
end

nDATA_COLO858 = DATA_COLO858./mean(DATA_COLO858(1,:,:),3);

% DATA_COLO858(Time<24,:,:) = [];
% Time(Time<24) = [];

DATA_mean = mean(DATA_COLO858,3);
DATA_std = std(DATA_COLO858,0,3);

%% Plot data



figure


for i = 1:6
    hold on
    errorbar(Time,DATA_mean(:,i),0.3*DATA_std(:,i),'LineWidth',2);
    % plot(Time,DATA_mean(:,i),'LineWidth',2)
    hold off
end

ax = gca;
ax.FontSize = 30;
ax.FontWeight = 'bold';
xlabel('Time (hours)')
ylabel('Cell count')

legend({'0\mu M','0.032 \mu M','0.1\mu M','0.32\mu M','1\mu M','3.2\mu M'},'Location','best')