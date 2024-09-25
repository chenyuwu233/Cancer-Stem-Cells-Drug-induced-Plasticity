%% This is a function to translate the parameters for 2-type hierarchy model 

function ret = hierarchy_translate(theta,max_cap)
    if max_cap == 0
        ret = theta;
    else
        c = theta(end);
        theta(end) = [];
        THETA = reshape(theta,[],2)';


        Drug_mat = THETA(1,4:end);
        Drug_mat = [Drug_mat;repmat(THETA(2,4:end),max_cap+1,1)];
        Cell_mat = [THETA(1,1:3),zeros(1,max_cap)];
        modified_Cell_vec = [0,THETA(2,2:3)];
        modified_Cell_mat = repmat(modified_Cell_vec,max_cap+1,1);
        modified_Cell_mat = [modified_Cell_mat,THETA(2,1)*eye(max_cap+1)];
        modified_Cell_mat(:,end) = [];
        Cell_mat = [Cell_mat;modified_Cell_mat];

        THETA = [Cell_mat,Drug_mat];
        ret = [reshape(THETA',1,[]),c];
    end

end