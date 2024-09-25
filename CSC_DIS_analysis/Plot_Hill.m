Conc = logspace(-3,1);

b = 0.85;

E = 0.125;

m = 2;

Hill =  b + (1-b)./(1+(Conc./E).^m);

semilogx(Conc,Hill,'LineWidth',3)
ax = gca;
ax.FontSize = 16;
% ax.FontWeight = 'bold';
xlabel('Drug concentration levels')
ylabel('Hill equation')
% title('Hill equation when $b >1$','Interpreter','latex')
