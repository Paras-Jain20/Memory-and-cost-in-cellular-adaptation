function [centered_moments, cond_prob_Sn_Sn_1,x_1_ss,x_0_ss,x_Sn_Sn_1, P_s] = marginal_resi_time_and_cond_prob_cell_state_stoc_env(N,q01,q10,P,varargin)
%% Stationary sitribution calculations
x = zeros(1,2*(N+1));
x(1) = 1;

x = x*P^1000;

%% Partitioning the State Transition Matrix for residence time calculations

Pb_01 = P(1:2*(floor(q01*N)+1),1:2*(floor(q01*N)+1));
Pb_bc_01 = P(1:2*(floor(q01*N)+1),2*(floor(q01*N)+1)+1:end);
Pbc_b_01 = P(2*(floor(q01*N)+1)+1:end,1:2*(floor(q01*N)+1));
Pbc_01 = P(2*(floor(q01*N)+1)+1:end,2*(floor(q01*N)+1)+1:end);

Pb_10 = P(1:2*(floor(q10*N)+1),1:2*(floor(q10*N)+1));
Pb_bc_10 = P(1:2*(floor(q10*N)+1),2*(floor(q10*N)+1)+1:end);
Pbc_b_10 = P(2*(floor(q10*N)+1)+1:end,1:2*(floor(q10*N)+1));
Pbc_10 = P(2*(floor(q10*N)+1)+1:end,2*(floor(q10*N)+1)+1:end);

%% Calculating P(Sn = 1) && P(Sn = 1 | Sn-1 = 1)
% Evolution of cell state with time assuming the starting cell state is 0

x_0 = zeros(2,size(Pb_01,1));
x_1 = zeros(2,size(Pbc_10,1));

xb = x(1:2*(floor(q01*N)+1));
x_0(1,:) = xb/sum(xb);

while(true)
    temp = x_1(1,:)*Pbc_10;
    x_1(2,:) = [temp(1:2*(floor(q01*N)+1)-size(Pb_10,1))  x_0(1,:)*Pb_bc_01+temp(2*(floor(q01*N)+1)+1-size(Pb_10,1):end)];

    temp = x_0(1,:)*Pb_01;
    x_0(2,:) = [temp(1:2*(floor(N*q10)+1))+x_1(1,:)*Pbc_b_10  temp(2*(floor(N*q10)+1)+1:end)];

    if(sum(abs(x_0(2,:)-x_0(1,:)))< 10^-9 && sum(abs(x_1(2,:)-x_1(1,:)))< 10^-9)
        x_1(1,:) = x_1(2,:);
        x_0(1,:) = x_0(2,:);
        break
    end

    x_0(1,:) = x_0(2,:);
    x_1(1,:) = x_1(2,:);
end

x_0_ss = x_0(1,:);
x_1_ss = x_1(1,:);

if(isempty(varargin))
    % for conditional prob assuming the cell state in 1
    x_0 = zeros(1,size(Pb_01,1));
    x_1(1,:) = x_1(1,:)/sum(x_1(1,:));

    temp = x_1(1,:)*Pbc_10;
    x_1(2,:) = [temp(1:2*(floor(q01*N)+1)-size(Pb_10,1))  x_0(1,:)*Pb_bc_01+temp(2*(floor(q01*N)+1)+1-size(Pb_10,1):end)];

    % for conditional prob assuming the cell state in 1
    if(sum(x_1_ss) > 0)
        cond_prob_Sn_Sn_1 = sum(x_1(2,:));

        x_0_old = zeros(1,size(Pb_01,1));
        x_1_old = x_1_ss/sum(x_1_ss);

        temp = x_1_old*Pbc_10;
        x_1_new = [temp(1:2*(floor(q01*N)+1)+1 -1-size(Pb_10,1)) x_0_old*Pb_bc_01+temp(2*(floor(q01*N)+1)+1-size(Pb_10,1):end)];
        % x_1_new = x_1_new/sum(x_1_new);
        x_Sn_Sn_1 = x_1_new;
    else
        cond_prob_Sn_Sn_1 = 0;
        x_Sn_Sn_1 = zeros(size(x_1_ss));
    end


    % determining transition probabilities between phenotypic states

    %P(S_n+1 = S_h|S_n = S_l) = P(k_n+1 = floor(q01*N)+1 |k_n = floor(q01*N), S_n = S_l) P(k_n = floor(q01*N)|S_n = S_l)
    %P(S_n+1 = S_l|S_n = S_h) = P(k_n+1 = floor(q10*N) |k_n = floor(q10*N)+1, S_n = S_h) P(k_n = floor(q10*N)+1|S_n = S_h)

    t_0_1 = P(2*(floor(q01*N)+1)-1,2*(floor(q01*N)+1)+2)*x_0_ss(end-1)/sum(x_0_ss) + P(2*(floor(q01*N)+1),2*(floor(q01*N)+1)+2)*x_0_ss(end)/sum(x_0_ss);
    t_0_0 = 1-t_0_1;
    t_1_0 = (P(2*(floor(q10*N)+1)+1,2*(floor(q10*N)+1)-1)*x_1_ss(1)/sum(x_1_ss) + P(2*(floor(q10*N)+1)+2,2*(floor(q10*N)+1)-1)*x_1_ss(2)/sum(x_1_ss) );
    t_1_1 = 1-t_1_0;

    P_s = [t_0_0 t_0_1; t_1_0 t_1_1];

    %% Residence time moment calculation
    centered_moments = zeros(4,2);
    uncentered_moments = zeros(4,2);

    Pb = P_s(1,1);
    Pb_bc = P_s(1,2);
    Pbc_b = P_s(2,1);
    Pbc = P_s(2,2);

    %%%%%%%%%%%%%%%%%% centered moments residence time in state 0 %%%%%%%%%%%%%%%%%%%%%%%%%

    v_0 = 1;
    I_Pb_inv = (1-Pb)^-1;

    % factorial moments
    fac_mom = zeros(4,1); % calculating factorial moment upto forth order
    temp = 1;
    temp1 = 1;

    for k = 1:4
        temp = I_Pb_inv*temp;
        if(k>1)
            temp1 = Pb*temp1;
        end
        fac_mom(k) = factorial(k)*v_0*temp1*temp;
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


    %%%%%%%%%%%%%%%%%% centered moments residence time in state 1 %%%%%%%%%%%%%%%%%%%%%%%%

    v_1 = 1;

    I_Pbc_inv = (1-Pbc)^-1;

    % factorial moments
    fac_mom = zeros(4,1); % calculating factorial moment upto forth order
    temp = 1;
    temp1 = 1;

    for k = 1:4
        temp = I_Pbc_inv*temp;
        if(k>1)
            temp1 = Pbc*temp1;
        end
        fac_mom(k) = factorial(k)*v_1*temp1*temp;
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
else
    cond_prob_Sn_Sn_1 = 0;
    x_Sn_Sn_1 = zeros(size(x_1_ss));
    P_s = zeros(2,2);
end
end