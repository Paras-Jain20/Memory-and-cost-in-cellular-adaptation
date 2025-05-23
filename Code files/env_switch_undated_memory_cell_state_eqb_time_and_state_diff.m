
function [adapt_time,state_difference] = env_switch_undated_memory_cell_state_eqb_time_and_state_diff(N,q01,q10,p1,p2)
% N = 10;
% q01 = 0.45;
% q10 = 0.15;
% p1 = 0;
% p2 = 0.4;

P1 = transition_matrix_P_cal(N,p1); 
P = transition_matrix_P_cal(N,p2);

% x = zeros(1,N+1);
% x(1) = 1;
% for i = 1:10000
% x = x*P1;
% end

[x_0_ss,x_1_ss,~] = undated_memory_cond_prob_cell_state(N,q01,q10,P1);

Pb_01 = P(1:floor(q01*N)+1,1:floor(q01*N)+1);
Pb_bc_01 = P(1:floor(q01*N)+1,floor(q01*N)+2:end);
Pbc_b_01 = P(floor(q01*N)+2:end,1:floor(q01*N)+1);
Pbc_01 = P(floor(q01*N)+2:end,floor(q01*N)+2:end);

Pb_10 = P(1:floor(q10*N)+1,1:floor(q10*N)+1);
Pb_bc_10 = P(1:floor(q10*N)+1,floor(q10*N)+2:end);
Pbc_b_10 = P(floor(q10*N)+2:end,1:floor(q10*N)+1);
Pbc_10 = P(floor(q10*N)+2:end,floor(q10*N)+2:end);

%% Calculating P(Sn = 1) && P(Sn = 1 | Sn-1 = 1)
% Evolution of cell state with time assuming the starting cell state is 0

x_0 = zeros(2,size(Pb_01,1));
x_1 = zeros(2,size(Pbc_10,1));

x_0(1,:) = x_0_ss;
x_1(1,:) = x_1_ss;


while(true)
    temp = x_1(1,:)*Pbc_10;
    x_1(2,:) = [temp(1:floor(q01*N)+1-size(Pb_10,1))  x_0(1,:)*Pb_bc_01+temp(floor(q01*N)+2-size(Pb_10,1):end)];

    temp = x_0(1,:)*Pb_01;
    x_0(2,:) = [temp(1:floor(N*q10)+1)+x_1(1,:)*Pbc_b_10  temp(floor(N*q10)+2:end)];

    if(sum(abs(x_0(2,:)-x_0(1,:)))< 10^-10 && sum(abs(x_1(2,:)-x_1(1,:)))< 10^-10)
        x_1(1,:) = x_1(2,:);
        x_0(1,:) = x_0(2,:);
        break
    end

    x_0(1,:) = x_0(2,:);
    x_1(1,:) = x_1(2,:);
%     adapt_time = adapt_time + 1;
end

x_0_ss_2 = x_0(1,:);
x_1_ss_2 = x_1(1,:);

prob_state_1_ss = sum(x_1_ss_2);
% state_difference = sum(x_1_ss_2) - sum(x_1_ss); 

x_0(1,:) = x_0_ss;
x_1(1,:) = x_1_ss;
adapt_time = 0;

prob_state_1 = zeros(1,2);
prob_state_1(1) = sum(x_1_ss);


while(true)
    temp = x_1(1,:)*Pbc_10;
    x_1(2,:) = [temp(1:floor(q01*N)+1-size(Pb_10,1))  x_0(1,:)*Pb_bc_01+temp(floor(q01*N)+2-size(Pb_10,1):end)];

    temp = x_0(1,:)*Pb_01;
    x_0(2,:) = [temp(1:floor(N*q10)+1)+x_1(1,:)*Pbc_b_10  temp(floor(N*q10)+2:end)];

    prob_state_1(2) = sum(x_1(2,:));
    
    if(abs(prob_state_1_ss - prob_state_1(2)) < 10^-10)
        break
    end
    x_0(1,:) = x_0(2,:);
    x_1(1,:) = x_1(2,:);
    prob_state_1(1) = prob_state_1(2);
    adapt_time = adapt_time + 1;
end

state_difference = prob_state_1(2)- sum(x_1_ss);
end

