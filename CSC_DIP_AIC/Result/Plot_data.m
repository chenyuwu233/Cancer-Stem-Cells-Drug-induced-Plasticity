%%

hold on

for i = 1:NC
    temp = squeeze(DATA(:,i,:))';
    errorbar(Time,mean(temp),std(temp),'LineWidth',3)

end

xlabel('Time')
ylabel('Total cell count')
legend(string(Conc),'Location','northwest')
title('without drug-induced plasticity')
ax = gca;
ax.FontSize = 23;
ax.FontWeight = 'bold';