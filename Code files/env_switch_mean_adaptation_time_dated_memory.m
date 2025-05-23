
function [centered_moments,prob_resi_time_state_0,prob_resi_time_state_1,prob_state_1_start_0,prob_state_1_start_1,prob_Sn_Sn_1_start_0,prob_Sn_Sn_1_start_1] = env_switch_mean_adaptation_time_dated_memory(N,q01,q10,p1,p2,total_time_steps)
% total_time_steps = 50;
centered_moments = zeros(4,2);
uncentered_moments = zeros(4,2); 

states = ff2n(N);
[states_sum,order_indx] = sort(sum(states,2),1,"ascend");
reordered_states = states(order_indx,:);

% % defining probability transition matrix P1
% k = 1;
% row =zeros(2*2^N,1);
% col = zeros(2*2^N,1);
% value = zeros(2*2^N,1);
% for row_indx = 1:2^N
% 
%     focal_state = reordered_states(row_indx,:);
%     updated_states = [0 focal_state(1:N-1); 1 focal_state(1:N-1)];
% 
%     col_indx = find(sum(reordered_states == updated_states(1,:),2) == N,1);
%     if(~isempty(col_indx))
%         row(k) = row_indx;
%         col(k) = col_indx;
%         value(k) = 1-p1;
%         k = k+1;
%     end
% 
%     col_indx = find(sum(reordered_states == updated_states(2,:),2) == N,1);
%     if(~isempty(col_indx))
%         row(k) = row_indx;
%         col(k) = col_indx;
%         value(k) = p1;
%         k = k+1;
%     end
% end
% 
% P1 = sparse(row,col,value);
% 
% x = zeros(1,2^N);
% x(1) = 1;
% for i = 1:10000
% x = x*P1;
% end

[x_0_ss,x_1_ss, ~] = dated_memory_cond_prob_cell_state(N,p1,q01,q10);

% defining probability transition matrix P
k = 1;
row =zeros(2*2^N,1);
col = zeros(2*2^N,1);
value = zeros(2*2^N,1);
for row_indx = 1:2^N

    focal_state = reordered_states(row_indx,:);
    updated_states = [0 focal_state(1:N-1); 1 focal_state(1:N-1)];

    col_indx = find(sum(reordered_states == updated_states(1,:),2) == N,1);
    if(~isempty(col_indx))
        row(k) = row_indx;
        col(k) = col_indx;
        value(k) = 1-p2;
        k = k+1;
    end
 
    col_indx = find(sum(reordered_states == updated_states(2,:),2) == N,1);
    if(~isempty(col_indx))
        row(k) = row_indx;
        col(k) = col_indx;
        value(k) = p2;
        k = k+1;
    end
end

P = sparse(row,col,value);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch_indx_01 = find(states_sum>(q01*N),1);
switch_indx_10 = find(states_sum>(q10*N),1);

%%%%%%%%%%%%%%%%%


Pb_01 = P(1:switch_indx_01-1,1:switch_indx_01-1);
Pb_bc_01 = P(1:switch_indx_01-1,switch_indx_01:end);
Pbc_b_01 = P(switch_indx_01:end,1:switch_indx_01-1);
Pbc_01 = P(switch_indx_01:end,switch_indx_01:end);

Pb_10 = P(1:switch_indx_10-1,1:switch_indx_10-1);
Pb_bc_10 = P(1:switch_indx_10-1,switch_indx_10:end);
Pbc_b_10 = P(switch_indx_10:end,1:switch_indx_10-1);
Pbc_10 = P(switch_indx_10:end,switch_indx_10:end);

%% Evolution of cell state with time assuming the starting cell state is 0
xb = x_0_ss; %x(1:switch_indx_01-1);
% xbc = x(floor(q01*N)+2:end);

x_0 = zeros(total_time_steps+1,size(Pb_01,1));
x_1 = zeros(total_time_steps,size(Pbc_10,1));

x_0(1,:) = xb/sum(xb);
% x_1(1,floor(q01*N)+2-size(Pb_10,1):end) = xbc;

prob_state_1_start_0 = zeros(1, total_time_steps);
prob_state_0_start_0 = zeros(1, total_time_steps);

prob_Sn_Sn_1_start_0 = zeros(1, total_time_steps);

for k = 1:total_time_steps
%     prob_resi_time_state_0(k) = v_0(k,:)*(prod2(:,:,k))*(eye(size(Pb_01))-Pb_01)*ones(size(Pb_01,1),1);
    prob_state_1_start_0(k) = sum(x_1(k,:));
    prob_state_0_start_0(k) = sum(x_0(k,:));


    temp = x_1(k,:)*Pbc_10; 
    x_1(k+1,:) = [temp(1:switch_indx_01-1-size(Pb_10,1)) x_0(k,:)*Pb_bc_01+temp(switch_indx_01-size(Pb_10,1):end)];

    temp = x_0(k,:)*Pb_01;
    x_0(k+1,:) = [temp(1:switch_indx_10-1)+x_1(k,:)*Pbc_b_10 temp(switch_indx_10:end)];

    % for conditional prob assuming the cell state in 1
    if(sum(x_1(k,:)) ~= 0)
        temp_x_0 = zeros(1,size(Pb_01,1));
        temp_x_1 = x_1(k,:)/sum(x_1(k,:));

        temp = temp_x_1*Pbc_10;
        temp_x_1 = [temp(1:switch_indx_01-1-size(Pb_10,1))  temp_x_0*Pb_bc_01+temp(switch_indx_01-size(Pb_10,1):end)];

        prob_Sn_Sn_1_start_0(k+1) = sum(temp_x_1)*prob_state_1_start_0(k);
    end

end

%% Evolution of cell state with time assuming the starting cell state is 1
% xb = x(1:floor(q10*N)+1);
xbc = x_1_ss; %x(switch_indx_10:end);

x_0 = zeros(total_time_steps+1,size(Pb_01,1));
x_1 = zeros(total_time_steps,size(Pbc_10,1));

% x_0(1,1:floor(N*q10)+1) = xb;
x_1(1,:) = xbc/sum(xbc);

prob_state_1_start_1 = zeros(1, total_time_steps);
prob_state_0_start_1 = zeros(1, total_time_steps);

prob_Sn_Sn_1_start_1 = zeros(1, total_time_steps);

for k = 1:total_time_steps
%     prob_resi_time_state_0(k) = v_0(k,:)*(prod2(:,:,k))*(eye(size(Pb_01))-Pb_01)*ones(size(Pb_01,1),1);
    prob_state_1_start_1(k) = sum(x_1(k,:));
    prob_state_0_start_1(k) = sum(x_0(k,:));

    temp = x_1(k,:)*Pbc_10; 
    x_1(k+1,:) = [temp(1:switch_indx_01-1-size(Pb_10,1)) x_0(k,:)*Pb_bc_01+temp(switch_indx_01-size(Pb_10,1):end)];

    temp = x_0(k,:)*Pb_01;
    x_0(k+1,:) = [temp(1:switch_indx_10-1)+x_1(k,:)*Pbc_b_10 temp(switch_indx_10:end)];

    % for conditional prob assuming the cell state in 1
    if(sum(x_1(k,:)) ~= 0)
        temp_x_0 = zeros(1,size(Pb_01,1));
        temp_x_1 = x_1(k,:)/sum(x_1(k,:));

        temp = temp_x_1*Pbc_10;
        temp_x_1 = [temp(1:switch_indx_01-1-size(Pb_10,1))  temp_x_0*Pb_bc_01+temp(switch_indx_01-size(Pb_10,1):end)];

        prob_Sn_Sn_1_start_1(k+1) = sum(temp_x_1)*prob_state_1_start_1(k);
    end

end

%% Probability of sojourn time in cell state 0 before switching to cell
% state 1
xb = x_0_ss; %x(1:switch_indx_01-1);
% xbc = x(switch_indx_01:end);

v_0 = xb/(sum(xb));
v_1 = zeros(total_time_steps,size(Pbc_10,1));

prob_resi_time_state_0 = zeros(1, total_time_steps);
joint_prob_residence_time_state_0_and_1 = zeros(total_time_steps, total_time_steps);

prod2 = eye(size(Pb_01)); 
for k = 1:total_time_steps
    prob_resi_time_state_0(k) = v_0*(prod2)*(eye(size(Pb_01))-Pb_01)*ones(size(Pb_01,1),1);
    v_1(k,switch_indx_01-size(Pb_10,1):end) = v_0*prod2*Pb_bc_01;


    prod1 = (eye(size(Pbc_10)))*(eye(size(Pbc_10))-Pbc_10)*ones(size(Pbc_10,1),1);
    for j = 1:total_time_steps
        joint_prob_residence_time_state_0_and_1(k,j) = v_1(k,:)*prod1;
        prod1 = Pbc_10*prod1;
    end
    prod2 = Pb_01*prod2;
end

% var_resi_time_state_0  = exp_resi_time_state_0_square - (mean_resi_time_state_0)^2;

%%% calculation of the higher order moments (particularly variance) of the
%%% residence time distribution

I_Pb_01_inv = sparse(((eye(size(Pb_01))-Pb_01)^-1));

% factorial moments
fac_mom = zeros(4,1); % calculating factorial moment upto forth order
temp = eye(size(I_Pb_01_inv));
temp1 = eye(size(Pb_01));

for k = 1:4
    temp = I_Pb_01_inv*temp;
    if(k>1)
        temp1 = Pb_01*temp1;
    end
    % fac_mom(k) = factorial(k)*v(end,:)*Pb^(k-1)*(eye(size(Pb))-Pb)^(-k)*eye(size(Pb,1),1);
    fac_mom(k) = factorial(k)*v_0*temp1*temp*ones(size(Pb_01,1),1);
end

for k = 1:size(uncentered_moments,1)
    for j = 1:k
        uncentered_moments(k,1) = uncentered_moments(k,1) + sterling_number(k,j)*fac_mom(j);
    end
end

% calculation of mean (mu) E(X) 
centered_moments(1,1) = uncentered_moments(1,1); 
% calculation of variance (sigma^2) E((X-mu)^2)
centered_moments(2,1) = uncentered_moments(2,1)-(uncentered_moments(1,1))^2;
% calculation of skewness E((X-mu)^3)/(sigma^3)
centered_moments(3,1)= (uncentered_moments(3,1) - 3*uncentered_moments(2,1)*uncentered_moments(1,1) + 2*uncentered_moments(1,1)^3)/(centered_moments(2,1)*sqrt(centered_moments(2,1)));
% calculation of kurtosis E((X-mu)^4)/(sigma^4)
centered_moments(4,1)= (uncentered_moments(4,1) - 4*uncentered_moments(3,1)*uncentered_moments(1,1) + 12*uncentered_moments(2,1)*uncentered_moments(1,1)^2 - 3*uncentered_moments(1,1)^4)/(centered_moments(2,1)^2);


%% Probability of sojourn time in cell state 1 before switching to cell
% state 0
% xb = x(1:switch_indx_10-1);
xbc = x_1_ss; %x(switch_indx_10:end);

v_1 = xbc/(sum(xbc));

prob_resi_time_state_1 = zeros(1, total_time_steps);

prod2 = eye(size(Pbc_10)); 
for k = 1:total_time_steps
    prob_resi_time_state_1(k) = v_1*(prod2)*(eye(size(Pbc_10))-Pbc_10)*ones(size(Pbc_10,1),1);
    prod2 = Pbc_10*prod2;
end

%%%%%%%%%%%%%%%%%% centered moments residence time in state 1 %%%%%%%%%%%%%%%%%%%%%%%%%

I_Pbc_10_inv = sparse(((eye(size(Pbc_10))-Pbc_10)^-1));

% factorial moments
fac_mom = zeros(4,1); % calculating factorial moment upto forth order
temp = eye(size(I_Pbc_10_inv));
temp1 = eye(size(Pbc_10));

for k = 1:4
    temp = I_Pbc_10_inv*temp;
    if(k>1)
        temp1 = Pbc_10*temp1;
    end
    % fac_mom(k) = factorial(k)*v(end,:)*Pb^(k-1)*(eye(size(Pb))-Pb)^(-k)*eye(size(Pb,1),1);
    fac_mom(k) = factorial(k)*v_1*temp1*temp*ones(size(Pbc_10,1),1);
end

for k = 1:size(uncentered_moments,1)
    for j = 1:k
        uncentered_moments(k,2) = uncentered_moments(k,2) + sterling_number(k,j)*fac_mom(j);
    end
end

% calculation of mean (mu) E(X) 
centered_moments(1,2) = uncentered_moments(1,2); 
% calculation of variance (sigma^2) E((X-mu)^2)
centered_moments(2,2) = uncentered_moments(2,2)-(uncentered_moments(1,2))^2;
% calculation of skewness E((X-mu)^3)/(sigma^3)
centered_moments(3,2)= (uncentered_moments(3,2) - 3*uncentered_moments(2,2)*uncentered_moments(1,2) + 2*uncentered_moments(1,2)^3)/(centered_moments(2,2)*sqrt(centered_moments(2,2)));
% calculation of kurtosis E((X-mu)^4)/(sigma^4)
centered_moments(4,2)= (uncentered_moments(4,2) - 4*uncentered_moments(3,2)*uncentered_moments(1,2) + 12*uncentered_moments(2,2)*uncentered_moments(1,2)^2 - 3*uncentered_moments(1,2)^4)/(centered_moments(2,2)^2);
end

% subplot(2,1,1)
% bar(prob_resi_time_state_0(end,:))
% xlim([0 50])
% 
% subplot(2,1,2)
% bar(prob_resi_time_state_1(end,:))
% xlim([0 50])

% here the v1 should be asymptotic prob of starting sojorn in state 1


