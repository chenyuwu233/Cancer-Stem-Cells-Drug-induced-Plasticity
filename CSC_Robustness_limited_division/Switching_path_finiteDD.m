%% Description:
% This is a path generating function for the switching model special for Differentiated cell lineage.
%
% Inputs (suppose we have 2 sub-types)
% init: 1 x 2 vector that includes initial cell numbers
% theta: 2 x 3 matrix that includes the birth rate, death rate, and
% switching rate among every sub-types
% Time: 1 x t vector that includes all the time points we collect the data.
%
% Output:
% ret: 2 x t vector that records the number of each type cells at every
% time points



function ret = Switching_path_finiteDD(init,theta,Time)
    SC = init(1); % Stem cell count
    DC = ones(1,init(2)); % Differentiated cell count vector
    DC_indi = ones(1,init(2)); % Indicator of DC resources. 
    DC_cap = 3;
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
        dc_rate   = DC_b*sum(DC(DC_indi<DC_cap))+(DC_d+DC_s)*sum(DC);
        rate      = sc_rate+dc_rate;
        t    = t - (1/rate) * log(rand);
        id = id + 1;
        if t > Time(nt)
            ret(:,nt) = [SC;sum(DC)];
            nt      = nt+1;
            if nt > length(Time)
                break
            end
        end
        if sum(DC)+SC > 1e8
            ret(:,nt:end) = [SC;sum(DC)];
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
                DC_indi = [DC_indi,1];
            end
        else % Differentiated Cell event
            dc_event = rand;
            dc_birth  = [0,DC_b*sum(DC(DC_indi<DC_cap))/dc_rate];
            dc_death  = [DC_b*sum(DC(DC_indi<DC_cap))/dc_rate,(1-DC_s*sum(DC)/dc_rate)];
            dc_switch = [(1-DC_s*sum(DC)/dc_rate),1];
            if dc_event > dc_birth(2)
                if dc_event>dc_switch(1) % Assume: every cell can de-differentiate
                    SC = SC+1;
                else % Assume: every cell can die
                    d_indi = search_grid(dc_event,dc_death,DC);
                    DC(d_indi) = DC(d_indi)-1;
                end
            else % Assume each lineage can only get DC_cap birth
                b_indi = search_grid(dc_event,dc_birth,DC(DC_indi<DC_cap));
                DC(b_indi) = DC(b_indi)+1;
                DC_indi(b_indi) = DC_indi(b_indi)+1;
            end
        end
    end
%     length(DC)
%     DC(DC_indi == DC_cap)  % Check whether we reach the DC_cap
end