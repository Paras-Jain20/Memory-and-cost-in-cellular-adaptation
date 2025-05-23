
function [adapt_time,state_difference] = env_switch_dated_memory_cell_state_eqb_time_and_state_diff(N,q01,q10,p1,p2)

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

% Evolution of cell state with time assuming the starting cell state is 0
% xb = x(1:switch_indx_01-1);
% xbc = x(floor(q01*N)+2:end);

x_0 = zeros(2,size(Pb_01,1));
x_1 = zeros(2,size(Pbc_10,1));


x_0(1,:) = x_0_ss;
x_1(1,:) = x_1_ss;


while(true)

    temp = x_1(1,:)*Pbc_10;
    x_1(2,:) = [temp(1:switch_indx_01-1-size(Pb_10,1)) x_0(1,:)*Pb_bc_01+temp(switch_indx_01-size(Pb_10,1):end)];

    temp = x_0(1,:)*Pb_01;
    x_0(2,:) = [temp(1:switch_indx_10-1)+x_1(1,:)*Pbc_b_10 temp(switch_indx_10:end)];

    if(sum(abs(x_0(2,:)-x_0(1,:)))< 10^-10 && sum(abs(x_1(2,:)-x_1(1,:)))< 10^-10)
        x_1(1,:) = x_1(2,:);
        x_0(1,:) = x_0(2,:);
        break
    end
    x_0(1,:) = x_0(2,:);
    x_1(1,:) = x_1(2,:);
%     adapt_time = adapt_time + 1;
end

prob_state_1_ss = sum(x_1(1,:));
% x_0_ss_2 = x_0(1,:);
% x_1_ss_2 = x_1(1,:);
% 
% state_difference = sum(x_1_ss_2) - sum(x_1_ss); 
prob_state_1 = zeros(1,2);
prob_state_1(1) = sum(x_1_ss);

x_0 = zeros(2,size(Pb_01,1));
x_1 = zeros(2,size(Pbc_10,1));
x_0(1,:) = x_0_ss;
x_1(1,:) = x_1_ss;

adapt_time = 0;
while(true)
    
    temp = x_1(1,:)*Pbc_10;
    x_1(2,:) = [temp(1:switch_indx_01-1-size(Pb_10,1)) x_0(1,:)*Pb_bc_01+temp(switch_indx_01-size(Pb_10,1):end)];

    temp = x_0(1,:)*Pb_01;
    x_0(2,:) = [temp(1:switch_indx_10-1)+x_1(1,:)*Pbc_b_10 temp(switch_indx_10:end)];

    prob_state_1(2) = sum(x_1(2,:));
    
    if(abs(prob_state_1(2) - prob_state_1_ss) < 10^-10)
        break
    end
    x_0(1,:) = x_0(2,:);
    x_1(1,:) = x_1(2,:);
    prob_state_1(1) = prob_state_1(2);
    adapt_time = adapt_time + 1;
end

state_difference = prob_state_1(2) - sum(x_1_ss);
end



