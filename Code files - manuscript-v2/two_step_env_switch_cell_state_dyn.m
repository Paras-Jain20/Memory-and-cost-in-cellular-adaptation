
function [prob_state_1_start_0,prob_Sn_Sh_Sn_1_Sh_start_0,prob_Sn_Sl_Sn_1_Sl_start_0] = two_step_env_switch_cell_state_dyn(N,q01,q10,end_env_pair,intermediate_env,total_time_steps,env_switch_time)
P_end = transition_matrix_P_cal(N,end_env_pair(2));
P = transition_matrix_P_cal(N,intermediate_env);

Pb_01 = P(1:floor(q01*N)+1,1:floor(q01*N)+1);
Pb_bc_01 = P(1:floor(q01*N)+1,floor(q01*N)+2:end);
Pbc_b_01 = P(floor(q01*N)+2:end,1:floor(q01*N)+1);
Pbc_01 = P(floor(q01*N)+2:end,floor(q01*N)+2:end);

Pb_10 = P(1:floor(q10*N)+1,1:floor(q10*N)+1);
Pb_bc_10 = P(1:floor(q10*N)+1,floor(q10*N)+2:end);
Pbc_b_10 = P(floor(q10*N)+2:end,1:floor(q10*N)+1);
Pbc_10 = P(floor(q10*N)+2:end,floor(q10*N)+2:end);

%% Evolution of cell state and P(Sn and Sn-1) with time assuming the starting cell state is 0

x_0 = zeros(total_time_steps+1,size(Pb_01,1));
x_1 = zeros(total_time_steps,size(Pbc_10,1));

x_0(1,1) = 1;

prob_state_1_start_0 = zeros(1, total_time_steps);
prob_state_0_start_0 = zeros(1, total_time_steps);

prob_Sn_Sh_Sn_1_Sh_start_0 = zeros(1, total_time_steps+1);
prob_Sn_Sh_Sn_1_Sh_start_0(1) = 0;

prob_Sn_Sl_Sn_1_Sl_start_0 = zeros(1, total_time_steps+1);
prob_Sn_Sl_Sn_1_Sl_start_0(1) = 1;

% the environment switches at k = 2 tp p > 0

for k = 1:env_switch_time
%     prob_resi_time_state_0(k) = v_0(k,:)*(prod2(:,:,k))*(eye(size(Pb_01))-Pb_01)*ones(size(Pb_01,1),1);
    prob_state_1_start_0(k) = sum(x_1(k,:));
    prob_state_0_start_0(k) = sum(x_0(k,:));

    temp = x_1(k,:)*Pbc_10; 
    x_1(k+1,:) = [temp(1:floor(q01*N)+1-size(Pb_10,1)) x_0(k,:)*Pb_bc_01+temp(floor(q01*N)+2-size(Pb_10,1):end)];

    temp = x_0(k,:)*Pb_01;
    x_0(k+1,:) = [temp(1:floor(N*q10)+1)+x_1(k,:)*Pbc_b_10 temp(floor(N*q10)+2:end)];

    % for conditional prob assuming the cell state in 1
    if(sum(x_1(k,:)) ~= 0)
        temp_x_0 = zeros(1,size(Pb_01,1));
        temp_x_1 = x_1(k,:)/sum(x_1(k,:));

        temp = temp_x_1*Pbc_10;
        temp_x_1 = [temp(1:floor(q01*N)+1-size(Pb_10,1))  temp_x_0*Pb_bc_01+temp(floor(q01*N)+2-size(Pb_10,1):end)];

        prob_Sn_Sh_Sn_1_Sh_start_0(k+1) = sum(temp_x_1)*prob_state_1_start_0(k);
    end

    % for conditional prob assuming the cell state in 0
    if(sum(x_0(k,:)) ~= 0)
        temp_x_0 = x_0(k,:)/sum(x_0(k,:));
        temp_x_1 = zeros(1,size(Pbc_10,1));

        temp = temp_x_0*Pb_01;
        temp_x_0 = [temp(1:floor(N*q10)+1)+temp_x_1*Pbc_b_10  temp(floor(N*q10)+2:end)];

        prob_Sn_Sl_Sn_1_Sl_start_0(k+1) = sum(temp_x_0)*prob_state_0_start_0(k);
    end
end


Pb_01 = P_end(1:floor(q01*N)+1,1:floor(q01*N)+1);
Pb_bc_01 = P_end(1:floor(q01*N)+1,floor(q01*N)+2:end);
Pbc_b_01 = P_end(floor(q01*N)+2:end,1:floor(q01*N)+1);
Pbc_01 = P_end(floor(q01*N)+2:end,floor(q01*N)+2:end);

Pb_10 = P_end(1:floor(q10*N)+1,1:floor(q10*N)+1);
Pb_bc_10 = P_end(1:floor(q10*N)+1,floor(q10*N)+2:end);
Pbc_b_10 = P_end(floor(q10*N)+2:end,1:floor(q10*N)+1);
Pbc_10 = P_end(floor(q10*N)+2:end,floor(q10*N)+2:end);

for k = env_switch_time+1:total_time_steps
    prob_state_1_start_0(k) = sum(x_1(k,:));
    prob_state_0_start_0(k) = sum(x_0(k,:));

    temp = x_1(k,:)*Pbc_10; 
    x_1(k+1,:) = [temp(1:floor(q01*N)+1-size(Pb_10,1)) x_0(k,:)*Pb_bc_01+temp(floor(q01*N)+2-size(Pb_10,1):end)];

    temp = x_0(k,:)*Pb_01;
    x_0(k+1,:) = [temp(1:floor(N*q10)+1)+x_1(k,:)*Pbc_b_10 temp(floor(N*q10)+2:end)];

    % for conditional prob assuming the cell state in 1
    if(sum(x_1(k,:)) ~= 0)
        temp_x_0 = zeros(1,size(Pb_01,1));
        temp_x_1 = x_1(k,:)/sum(x_1(k,:));

        temp = temp_x_1*Pbc_10;
        temp_x_1 = [temp(1:floor(q01*N)+1-size(Pb_10,1))  temp_x_0*Pb_bc_01+temp(floor(q01*N)+2-size(Pb_10,1):end)];

        prob_Sn_Sh_Sn_1_Sh_start_0(k+1) = sum(temp_x_1)*prob_state_1_start_0(k);
    end

        % for conditional prob assuming the cell state in 0
    if(sum(x_0(k,:)) ~= 0)
        temp_x_0 = x_0(k,:)/sum(x_0(k,:));
        temp_x_1 = zeros(1,size(Pbc_10,1));

        temp = temp_x_0*Pb_01;
        temp_x_0 = [temp(1:floor(N*q10)+1)+temp_x_1*Pbc_b_10  temp(floor(N*q10)+2:end)];

        prob_Sn_Sl_Sn_1_Sl_start_0(k+1) = sum(temp_x_0)*prob_state_0_start_0(k);
    end
end

% accounting for the one time step delay in manifesting change in memory to
% change in cell state; Now the probability of cell state are aligned with
% change of the env
prob_state_1_start_0 = [0 prob_state_1_start_0];
prob_Sn_Sl_Sn_1_Sl_start_0 = [1 prob_Sn_Sl_Sn_1_Sl_start_0]; 
prob_Sn_Sh_Sn_1_Sh_start_0 = [0 prob_Sn_Sh_Sn_1_Sh_start_0];
end

