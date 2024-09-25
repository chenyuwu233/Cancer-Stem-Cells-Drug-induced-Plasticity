 
function ret = plot_sp(theta,Conc,cmd)
    ax = gca;
    True_sp = [];
    Theta   = coeff_reshape(theta(1:end-1),cmd,'Like',2);
    % Theta   = sortrows(Theta,5);
    for j = 2:length(Conc)
            A_j = Drug_A(Theta,Conc(j),cmd);
            sp  = get_stable_p(A_j);
        %     lam
        %     pause
            True_sp = [True_sp,sp(1)];
    end
    
    % Sp_prc = prctile(Sp_est,[2.5,97.5]);
    % Sp_prc_neg = median(Sp_est) - min(Sp_prc);
    % Sp_prc_pos = max(Sp_prc) - median(Sp_est);
    % errorbar([1,2,3,4,5,6,7,8,9],median(Sp_est),Sp_prc_neg,Sp_prc_pos)
    % boxplot(Sp_est,{'C1','C2','C3','C4','C5','C6','C7','C8','C9'})
    Conc_temp = Conc;
    Conc_temp(1) = [];
%     ax.XTickLabel(Conc_temp)
%     Conc(1) = 1e-6;
    ylim([0,1])
    xlabel('Concentration levels')
    ylabel('Stable proportion')
    hold on
    plot(Conc_temp,True_sp,'-o')
    ax.XScale = 'log';
    ax.XLim = [min(Conc_temp),max(Conc_temp)];
    ax.XTickLabel = Conc_temp;
%     xticklabels(Conc_temp)
    ret = True_sp;
end