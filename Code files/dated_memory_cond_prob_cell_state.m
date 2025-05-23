function [x_0_ss,x_1_ss, cond_prob_Sn_Sn_1] = dated_memory_cond_prob_cell_state(N,p,q01,q10)
%% State Transition Matrix
states = ff2n(N);
[states_sum,order_indx] = sort(sum(states,2),1,"ascend");
reordered_states = states(order_indx,:);

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
        value(k) = 1-p;
        k = k+1;
    end

    col_indx = find(sum(reordered_states == updated_states(2,:),2) == N,1);
    if(~isempty(col_indx))
        row(k) = row_indx;
        col(k) = col_indx;
        value(k) = p;
        k = k+1;
    end
end

P = sparse(row,col,value);

%% Partitioning the State Transition Matrix for residence time calculations

switch_indx_01 = find(states_sum>(q01*N),1);
switch_indx_10 = find(states_sum>(q10*N),1);

x = zeros(1,2^N);
x(1) = 1;
for i = 1:10000
    x = x*P;
end

Pb_01 = P(1:switch_indx_01-1,1:switch_indx_01-1);
Pb_bc_01 = P(1:switch_indx_01-1,switch_indx_01:end);
Pbc_b_01 = P(switch_indx_01:end,1:switch_indx_01-1);
Pbc_01 = P(switch_indx_01:end,switch_indx_01:end);

Pb_10 = P(1:switch_indx_10-1,1:switch_indx_10-1);
Pb_bc_10 = P(1:switch_indx_10-1,switch_indx_10:end);
Pbc_b_10 = P(switch_indx_10:end,1:switch_indx_10-1);
Pbc_10 = P(switch_indx_10:end,switch_indx_10:end);

xb = x(1:switch_indx_10-1);
xbc = x(switch_indx_10:end);

%% Calculating P(Sn = 1) && P(Sn = 1 | Sn-1 = 1)
if(q01 ~= q10)

    % Evolution of cell state with time assuming the starting cell state is 0
    % xb = x(1:switch_indx_01-1);
    % xbc = x(floor(q01*N)+2:end);

    x_0 = zeros(2,size(Pb_01,1));
    x_1 = zeros(2,size(Pbc_10,1));
    
    % x_0(1,1) = 1;
    if(p ~= 1)
        xb = x(1:switch_indx_01-1);
        x_0(1,:) = xb/sum(xb);
    else
        xbc = x(switch_indx_10:end);
        x_1(1,:) = xbc/sum(xbc);
    end

    % x_0(1,:) = xb/sum(xb);
    % x_1(1,floor(q01*N)+2-size(Pb_10,1):end) = xbc;
    count = 1;
    while(true)
        temp = x_1(1,:)*Pbc_10;
        x_1(2,:) = [temp(1:switch_indx_01-1-size(Pb_10,1)) x_0(1,:)*Pb_bc_01+temp(switch_indx_01-size(Pb_10,1):end)];

        temp = x_0(1,:)*Pb_01;
        x_0(2,:) = [temp(1:switch_indx_10-1)+x_1(1,:)*Pbc_b_10 temp(switch_indx_10:end)];

        if(sum(abs(x_0(2,:)-x_0(1,:)))< 10^-6 && sum(abs(x_1(2,:)-x_1(1,:)))< 10^-6)
            x_1(1,:) = x_1(2,:);
            x_0(1,:) = x_0(2,:);
            break
        end
        x_0(1,:) = x_0(2,:);
        x_1(1,:) = x_1(2,:);
        count = count + 1;
    end
    
    x_0_ss = x_0(1,:);
    x_1_ss = x_1(1,:);
    % for conditional prob assuming the cell state in 1
    x_0 = zeros(1,size(Pb_01,1));
    x_1(1,:) = x_1(1,:)/sum(x_1(1,:));

    temp = x_1(1,:)*Pbc_10;
    x_1(2,:) = [temp(1:switch_indx_01-1-size(Pb_10,1)) x_0(1,:)*Pb_bc_01+temp(switch_indx_01-size(Pb_10,1):end)];

    cond_prob_Sn_Sn_1 = sum(x_1(2,:));

else
    % xb = x(1:switch_indx_10-1);
    xbc = x(switch_indx_10:end);

    x_0_ss = x(1:switch_indx_01-1);
    x_1_ss = x(switch_indx_10:end);

    cond_prob_Sn_Sn_1 = (1-sum(xbc/(sum(xbc))*Pbc_b_10));
end
end