init   = 1e4;
s      = 2;
Conc = [0,4,8];
NT   = 5;
Time = [2,6,12,24,48];
NC   = length(Conc);
NR   = 3;
cmd  = 'CSC_DIS';
lam_s_range       = [0,0.1];
lam_d_range       = [0,0.5];


Cell_SC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figures 2A');
Cell_TC = readcell('CPX paper (Figures 1E, 1F, 2A and 2B data).xlsx','Sheet','Figure 2B');
DMSO_TC = readcell('CPX paper (Figures 1E and 2B data).xlsx','Sheet','Figure 2B');


DATA_SC = [cell2mat(Cell_SC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_SC(3:5,11:16))];
DATA_TC = [cell2mat(Cell_TC(3:5,2:9)),[nan;nan;nan],cell2mat(Cell_TC(3:5,11:16))];
DATA_DMSO = [cell2mat(DMSO_TC(9:11,2:9)),[nan;nan;nan],cell2mat(DMSO_TC(9:11,11:16))];
DATA_DMSO_avg = mean(DATA_DMSO,1);

DATA    = zeros(NT,NC,NR);
DATA_sc = zeros(NT,NC,NR);
for i = 1:5
    k = i;
    DATA(i,:,:) = DATA_DMSO_avg(3*k-2:3*k).*DATA_TC(:,3*k-2:3*k)/100;
    DATA_sc(i,:,:) = squeeze(DATA(i,:,:)).*DATA_SC(:,3*k-2:3*k)/100;
end

E_data = mean(DATA(:,1,:),3,'omitnan');

Fit_5 = fit(Time',E_data,'exp1');




