%%


% tic
% opt_fval = get_like(DATA,theta,Time,Conc,NR,NC,NT,s,cmd);
% t_like = toc



tic 
opt_fval2 = get_like_alt(DATA,theta,Time,Conc,NR,NC,NT,s,cmd);
t_like_alt = toc