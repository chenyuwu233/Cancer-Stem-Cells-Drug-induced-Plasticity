% Inputs (suppose we have 2 sub-types)
%  - init: 1 x 2 vector that includes initial cell numbers
%  - theta: 2 x 3 matrix that includes the birth rate, death rate, and
%           switching rate among every sub-types
%  - Time: 1 x t vector that includes all the time points we collect the data.
%  - max_cap: maximum generation capacity
%  - cmd: switching dynamic assumption: 1. switching, 2. asymmetric
%
% Output:
% ret: 2 x t vector that records the number of each type cells at every
% time points



function ret = Switching_path_finiteDD_alt(init,theta,Time,max_cap,cmd)
    SC = init(1); % Stem cell count
    DC = ones(1,init(2)); % Differentiated cell count vector, both live and dead cell.
    DC_cap = max_cap;
    DC_indi = zeros(2^DC_cap,init(2)); % Indicator of DC generation
    DC_indi(DC,:) = 1;
    if size(theta,1) ~= 2 
        size(theta,1)
        pause
    end
    SC_b = theta(1,1);
    SC_d = theta(1,2);
    SC_s = theta(1,3);
    SC_rate = sum(theta(1,:));
    DC_b = theta(2,1);
    DC_d = theta(2,2);
    DC_s = theta(2,3);
    t  = Time(1);
    ret = zeros(2,length(Time));
    ret(:,1) = init';
    nt = 2;
    id = 0;
    while t < Time(nt)&&SC+sum(DC)>0
        sc_rate   = SC_rate*SC;
        dc_rate   = DC_b*sum(DC_indi>0 & DC_indi<=DC_cap,'all')+(DC_d+DC_s)*sum(DC_indi>0,"all");
        rate      = sc_rate+dc_rate;
        t    = t - (1/rate) * log(rand);
        id = id + 1;
        if t > Time(nt)
            ret(:,nt) = [SC;sum(DC_indi>0,'all')];
%             DC_indi
%             pause
            nt      = nt+1;
            if nt > length(Time)
                break
            end
        end
        if sum(DC_indi>0,'all')+SC > 1e8
            ret(:,nt:end) = [SC;sum(DC_indi>0,'all')];
            break
        end
        if sum(DC_indi>0,'all')+SC == 0
            for lm = nt:length(Time)
                ret(:,lm) = [0;0];
            end
            break
        end

        event  = rand;
        if event < sc_rate/rate % Stem Cell event
            sc_event = rand;
            if sc_event < SC_b/SC_rate % Stem Cell birth
                SC = SC+1;
            elseif SC_b/SC_rate <= sc_event && sc_event < (SC_b+SC_d)/SC_rate % Stem Cell death
                SC = SC-1;
            else % Stem Cell differentiate
                DC = [DC,1];
                new_cell = zeros(2^DC_cap,1);
                new_cell(1) = 1;
                DC_indi = [DC_indi,new_cell];
            end
        else % Differentiated Cell event

            dc_event = rand;
            dc_birth  = [0,DC_b*sum(DC_indi>0 & DC_indi<=DC_cap,'all')/dc_rate];
            dc_death  = [DC_b*sum(DC_indi>0 & DC_indi<=DC_cap,'all')/dc_rate,(1-DC_s*sum(DC_indi>0,'all')/dc_rate)];
            dc_switch = [(1-DC_s*sum(DC_indi>0,'all')/dc_rate),1];
            if dc_event > dc_birth(2)
                if dc_event>dc_switch(1) % Assume: every cell can de-differentiate, de-differentiation does not affect the differentiated cell.
                    SC = SC+1;
                    if strcmp(cmd,'switching')
                    % Assume switching
                        if size(DC_indi,1) == 1
                            lineage_indi = get_rand_index(DC_indi);
                            cell_indi = 1;
                        else
                            lineage_indi = get_rand_index(sum(DC_indi>0));
                            cell_indi = get_rand_index(DC_indi(:,lineage_indi)>0); % Cell indi should be less than the first zero in DC_indi
                        end
                        DC_indi(cell_indi,lineage_indi) = -1;
                    end
                else % Assume: every cell can die
                    if size(DC_indi,1) == 1
                        lineage_indi = get_rand_index(DC_indi);
                        cell_indi = 1;
                    else
                        lineage_indi = get_rand_index(sum(DC_indi>0));
                        cell_indi = get_rand_index(DC_indi(:,lineage_indi)>0); % Cell indi should be less than the first zero in DC_indi
                    end
                    DC_indi(cell_indi,lineage_indi) = -1;
                end
            else % Assume each lineage can only get DC_cap birth
                
                if size(DC_indi,1) == 1
                    cell_indi = 1;
                    lineage_indi = get_rand_index(DC_indi);
                    DC(lineage_indi) = DC(lineage_indi)+1;
                else
                    try
                        lineage_indi = get_rand_index(sum(DC_indi>0&DC_indi<=DC_cap));
                    catch
                        sc_rate
                        rate
                        event
                        SC
                        DC
                        keyboard
                    end
                    DC(lineage_indi) = DC(lineage_indi)+1;
                    cell_indi = get_rand_index(DC_indi(:,lineage_indi)>0&DC_indi(:,lineage_indi)<=DC_cap);
                    if DC_indi(DC(lineage_indi),lineage_indi) ~=0
                        DC_indi(:,lineage_indi)
                        pause
                    end
                end
                DC_indi(cell_indi,lineage_indi) = DC_indi(cell_indi,lineage_indi)+1;
                DC_indi(DC(lineage_indi),lineage_indi) = DC_indi(cell_indi,lineage_indi); % Duplicate the cell and information.
            end
        end
    end
end