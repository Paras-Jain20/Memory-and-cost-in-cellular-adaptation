function [residence_time_0_state,residence_time_1_state] = residence_time_cal_cell_state_based(cell_state)
%%%%%%%%%%%%%%%%%%%%%%%  Residence time calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
residence_time_0_state= zeros(length(cell_state),1);
residence_time_1_state= zeros(length(cell_state),1);
resi_indx_0 = 1;
resi_indx_1 = 1;

if(cell_state(1) == 0)
    start_time_0 = 1;
elseif(cell_state(1) == 1)
    start_time_1 = 1;
end

for t = 1:length(cell_state)-1
    if(cell_state(t) == 1 && cell_state(t+1) == 0)
        start_time_0 = t+1;
        end_time_1 = t;
        residence_time_1_state(resi_indx_1) = end_time_1 - start_time_1 + 1;
        resi_indx_1 = resi_indx_1 + 1;
    elseif(cell_state(t) == 0 && cell_state(t+1) == 1)
        start_time_1 = t+1;
        end_time_0 = t;
        residence_time_0_state(resi_indx_0) = end_time_0 - start_time_0 + 1;
        resi_indx_0 = resi_indx_0 + 1;
    elseif(t+1 == length(cell_state))
        if(cell_state(t+1) == 0 && resi_indx_0 == 1)
            end_time_0 = t+1;
            residence_time_0_state(resi_indx_0) = end_time_0 - start_time_0 + 1;
        elseif(cell_state(t+1) == 1 && resi_indx_1 == 1)
            end_time_1 = t+1;
            residence_time_1_state(resi_indx_1) = end_time_1 - start_time_1 + 1;
        end
    end
end

residence_time_0_state = nonzeros(residence_time_0_state);
residence_time_1_state = nonzeros(residence_time_1_state);

end