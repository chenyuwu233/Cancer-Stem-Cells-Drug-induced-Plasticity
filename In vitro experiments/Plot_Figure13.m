load('Result/In_vitro_COLO858_60h_AIC.mat')


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


%% Compute the mean
theta = opt_xx_pe_DIP_Asy;
c = theta(end)';


theta(end) = [];
tv = theta(end-1:end); % Variable to control the delay effect 
theta(end-1:end) = [];


Theta = reshape(theta,[],s)';
p_init     = Theta(:,1)';
Theta(:,1) = [];
b    = Theta(:,1)';

PE_mean = zeros(NC,NT);
p = p_init;



for i = 1:NC
    init_M = squeeze(DATA(1,i,:))*p;

    init = mean(init_M);
    PE_mean(i,1) = sum(init);
    
    Mij = get_Mean_TD(Theta,Conc(i),tv,Time_ss-Time_ss(1));
    for j = 2:NT
        Meani = sum(init*Mij(:,:,j-1));
        PE_mean(i,j) = Meani;
    end

end





%% Plot estimation
Time_60_full = Time(find(Time == 2):find(Time == 60));

DATA_COLO858_60_full = DATA_COLO858(find(Time == 2):find(Time == 60),:,:);


h = tiledlayout(2,3);

for i = 1:6

    ax = nexttile;
    
    
    
    % pDATA_mean = nDATA_mean(1:idx,:);
    % pDATA_std = nDATA_std(1:idx,:);
    

    hold on
    

    errorbar(Time_60_full,mean(DATA_COLO858_60_full(:,i,:),3),std(DATA_COLO858_60_full(:,i,:),0,3));
    plot(Time_ss,PE_mean(i,:),'LineWidth',3);
    % plot(Time_ss,sen_pop(i,:),'LineWidth',3,'LineStyle','-.')
    % plot(Time_ss,res_pop(i,:),'LineWidth',3,'LineStyle','-.')
    ax.FontSize = 20;
    ax.FontWeight = 'bold';
    tc = strcat('Fitting plot under: ',num2str(Conc(i)),' \mu M of Vemurafenib');
    title(tc)
    xlabel('Time (hours)')
    ylabel('Cell count')
    legend({'Data','Fitted mean'},'Location','best')
    
end