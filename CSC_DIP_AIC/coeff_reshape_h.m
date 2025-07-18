%% Description:
%  This is a function that reshape the theta into the desired form
%  This version consider the Hill coefficient with the Multi-selective model.
%  
%  Input (suppose we have s sub-types):
%  - theta: 1 x d vector that record all the parameter about the
%  sub-population
%  - cmd1: model command
%  - cmd2: usage command
%  - s: number of sub-populations
%
%  Output:
%  - Theta: desired parameter format
%  - p: initial proportion of each sub-population.




function [Theta,p] = coeff_reshape_h(theta,cmd1,cmd2,s)
    switch cmd2
        case 'Gen'
            switch cmd1 
                case {'Linear_switching','Linear_switching_EP'}
                    if length(theta) ~= s*(s+3)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+3,s])';
                    p     = Theta(:,1)';
                    Theta(:,1) = [];
                case {'Single_Broad','Hill2_switching_2','Single_Broad_d','Hill2_death_2','Multi_Selective','Hill2_death_2_i',...
                        'Single_Broad_EP','Hill2_switching_2_EP','Hill2_death_EP','Hill2_death_2_EP','Multi_Selective_EP','Hill2_death_2_i_EP'}
                    if length(theta) ~= s*(s+5)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+5,s])';
                    p     = Theta(:,1)';
                    Theta(:,1) = [];
%                 case {}
%                     if length(theta) ~= s*(s+4)
%                         warning('The dimension is not match.')
%                         pause
%                     end
%                     Theta = reshape(theta,[s+4,s])';
%                     p     = Theta(:,1)';
%                     Theta(:,1) = [];
                case {'Hill2_switching_death','Hill2_switching_death_EP'}
                    if length(theta) ~= s*(s+6)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+6,s])';
                    p     = Theta(:,1)';
                    Theta(:,1) = [];
                case {'Hill2_switching_Hill2_death','Hill2_switching_Hill2_switching',...
                        'Hill2_death_Hill2_death','Hill2_switching_Hill2_death_EP','Hill2_switching_Hill2_switching_EP',...
                        'Hill2_death_Hill2_death_EP'}
                    if length(theta) ~= s*(s+6)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+6,s])';
                    p     = Theta(:,1)';
                    Theta(:,1) = [];
                case {'Hill2_switching_Hill2_death_i','Multi_Broad',...
                        'Hill2_death_Hill2_death_i','Hill2_switching_Hill2_death_i_EP','Multi_Broad_EP',...
                        'Hill2_death_Hill2_death_i_EP'}
                    if length(theta) ~= s*(s+6)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+6,s])';
                    p     = Theta(:,1)';
                    Theta(:,1) = [];
            end


        case 'Like'
            switch cmd1
                case {'Linear_switching','Linear_switching_EP'}
                    if length(theta) ~= s*(s+3)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+3,s])';
                    p     = Theta(:,1)';
                    Theta(:,1) = [];
                case {'Single_Broad','Hill2_switching_2','Single_Broad_d','Hill2_death_2','Multi_Selective','Hill2_death_2_i',...
                        'Single_Broad_EP','Hill2_switching_2_EP','Hill2_death_EP','Hill2_death_2_EP','Multi_Selective_EP','Hill2_death_2_i_EP'}
                    if length(theta) ~= s*(s+5)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+5,s])';
                    p     = Theta(:,1)';
                    Theta(:,1) = [];
%                 case {}
%                     if length(theta) ~= s*(s+4)
%                         warning('The dimension is not match.')
%                         pause
%                     end
%                     Theta = reshape(theta,[s+4,s])';
%                     p     = Theta(:,1)';
%                     Theta(:,1) = [];
                case {'Hill2_switching_death','Hill2_switching_death_EP'}
                    if length(theta) ~= s*(s+6)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+6,s])';
                    p     = Theta(:,1)';
                    Theta(:,1) = [];
                case {'Hill2_switching_Hill2_death','Hill2_switching_Hill2_switching',...
                        'Hill2_death_Hill2_death','Hill2_switching_Hill2_death_EP','Hill2_switching_Hill2_switching_EP',...
                        'Hill2_death_Hill2_death_EP'}
                    if length(theta) == s*(s+6)
                        Theta = reshape(theta,[s+6,s])';
                        p     = Theta(:,1)';
                        Theta(:,1) = [];
                    elseif length(theta) == s*(s+6)-1
                        p     = theta(1:s+6:end-s-6);
                        theta(1:s+6:end-s-6) = [];
                        Theta = reshape(theta,[s+5,s])';
                        p     = [p,1-sum(p)];
                    else
                        warning('The dimension is not match.')
                        pause
                    end
                case {'Hill2_switching_Hill2_death_i','Multi_Broad',...
                        'Hill2_death_Hill2_death_i','Hill2_switching_Hill2_death_i_EP','Multi_Broad_EP',...
                        'Hill2_death_Hill2_death_i_EP'}
                    if length(theta) == s*(s+6)
                        Theta = reshape(theta,[s+6,s])';
                        p     = Theta(:,1)';
                        Theta(:,1) = [];
                    elseif length(theta) == s*(s+6)-1
                        p     = theta(1:s+6:end-s-6);
                        theta(1:s+6:end-s-6) = [];
                        Theta = reshape(theta,[s+5,s])';
                        p     = [p,1-sum(p)];
                    else
                        warning('The dimension is not match.')
                        pause
                    end
            end
        case 'sp'
            p = [];
            switch cmd1 
                case {'Linear_switching','Linear_switching_EP'}
                    if length(theta) ~= s*(s+2)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+2,s])';
                case {'Single_Broad','Hill2_switching_2','Single_Broad_d','Hill2_death_2','Multi_Selective','Hill2_death_2_i',...
                        'Single_Broad_EP','Hill2_switching_2_EP','Hill2_death_EP','Hill2_death_2_EP','Multi_Selective_EP','Hill2_death_2_i_EP'}
                    if length(theta) ~= s*(s+4)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+4,s])';
%                 case {}
%                     if length(theta) ~= s*(s+3)
%                         warning('The dimension is not match.')
%                         pause
%                     end
%                     Theta = reshape(theta,[s+3,s])';
                case {'Hill2_switching_death','Hill2_switching_death_EP'}
                    if length(theta) ~= s*(s+5)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+5,s])';
                case {'Hill2_switching_Hill2_death','Hill2_switching_Hill2_switching',...
                        'Hill2_death_Hill2_death','Hill2_switching_Hill2_death_EP','Hill2_switching_Hill2_switching_EP',...
                        'Hill2_death_Hill2_death_EP'}
                    if length(theta) ~= s*(s+5)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+5,s])';
                case {'Hill2_switching_Hill2_death_i','Multi_Broad',...
                        'Hill2_death_Hill2_death_i','Hill2_switching_Hill2_death_i_EP','Multi_Broad_EP',...
                        'Hill2_death_Hill2_death_i_EP'}
                    if length(theta) ~= s*(s+5)
                        warning('The dimension is not match.')
                        pause
                    end
                    Theta = reshape(theta,[s+5,s])';
            end  
    end
end