function [x_0_ss,x_1_ss,cond_prob_Sn_Sn_1] = undated_memory_cond_prob_cell_state(N,q01,q10,P)
%% Stationary sitribution calculations
p = P(1,2);

x = zeros(1,N+1);
if(p == 1)
    x(end) = 1;
elseif(p == 0)
    x(1) = 1;
else
    for x_indx = 1:N+1
        x(x_indx) = nchoosek(N,x_indx-1)*(p)^(x_indx-1)* (1-p)^(N-(x_indx-1));
    end
end

%% Partitioning the State Transition Matrix for residence time calculations

Pb_01 = P(1:floor(q01*N)+1,1:floor(q01*N)+1);
Pb_bc_01 = P(1:floor(q01*N)+1,floor(q01*N)+2:end);
Pbc_b_01 = P(floor(q01*N)+2:end,1:floor(q01*N)+1);
Pbc_01 = P(floor(q01*N)+2:end,floor(q01*N)+2:end);

Pb_10 = P(1:floor(q10*N)+1,1:floor(q10*N)+1);
Pb_bc_10 = P(1:floor(q10*N)+1,floor(q10*N)+2:end);
Pbc_b_10 = P(floor(q10*N)+2:end,1:floor(q10*N)+1);
Pbc_10 = P(floor(q10*N)+2:end,floor(q10*N)+2:end);


%% Calculating P(Sn = 1) && P(Sn = 1 | Sn-1 = 1)
if(q01 ~= q10)
    % Evolution of cell state with time assuming the starting cell state is 0
    
    x_0 = zeros(2,size(Pb_01,1));
    x_1 = zeros(2,size(Pbc_10,1));

    % x_0(1,1) = 1;
    if(p ~= 1)
        xb = x(1:floor(q01*N)+1);
        x_0(1,:) = xb/sum(xb);
    else
        xbc = x(floor(q10*N)+2:end);
        x_1(1,:) = xbc/sum(xbc);
    end

    while(true)
        temp = x_1(1,:)*Pbc_10;
        x_1(2,:) = [temp(1:floor(q01*N)+1-size(Pb_10,1))  x_0(1,:)*Pb_bc_01+temp(floor(q01*N)+2-size(Pb_10,1):end)];

        temp = x_0(1,:)*Pb_01;
        x_0(2,:) = [temp(1:floor(N*q10)+1)+x_1(1,:)*Pbc_b_10  temp(floor(N*q10)+2:end)];

        if(sum(abs(x_0(2,:)-x_0(1,:)))< 10^-6 && sum(abs(x_1(2,:)-x_1(1,:)))< 10^-6)
            x_1(1,:) = x_1(2,:);
            x_0(1,:) = x_0(2,:);
            break
        end

        x_0(1,:) = x_0(2,:);
        x_1(1,:) = x_1(2,:);
    end
    
    x_0_ss = x_0(1,:);
    x_1_ss = x_1(1,:);
    
    % for conditional prob assuming the cell state in 1
    x_0 = zeros(1,size(Pb_01,1));
    x_1(1,:) = x_1(1,:)/sum(x_1(1,:));

    temp = x_1(1,:)*Pbc_10;
    x_1(2,:) = [temp(1:floor(q01*N)+1-size(Pb_10,1))  x_0(1,:)*Pb_bc_01+temp(floor(q01*N)+2-size(Pb_10,1):end)];

    cond_prob_Sn_Sn_1 = sum(x_1(2,:));

else

    xbc = x(floor(N*q10)+2:end);
    
    x_0_ss = x(1:floor(N*q01)+1);
    x_1_ss = x(floor(N*q10)+2:end);

    cond_prob_Sn_Sn_1 = (1-sum(xbc/(sum(xbc))*Pbc_b_10));
end
end