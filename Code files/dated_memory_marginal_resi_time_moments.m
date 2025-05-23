function [centered_moments] = dated_memory_marginal_resi_time_moments(N,p,q01,q10)
centered_moments = zeros(4,2);
uncentered_moments = zeros(4,2);

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
% Pbc_b_01 = P(switch_indx_01:end,1:switch_indx_01-1);
% Pbc_01 = P(switch_indx_01:end,switch_indx_01:end);

Pb_10 = P(1:switch_indx_10-1,1:switch_indx_10-1);
% Pb_bc_10 = P(1:switch_indx_10-1,switch_indx_10:end);
Pbc_b_10 = P(switch_indx_10:end,1:switch_indx_10-1);
Pbc_10 = P(switch_indx_10:end,switch_indx_10:end);

% xb = x(1:switch_indx_10-1);
xbc = x(switch_indx_10:end);

%% Residence time moments calculations

%%%%%%%%%%%%%%%%%% centered moments residence time in state 0 %%%%%%%%%%%%%%%%%%%%%%%%%

I_Pb_01_inv = sparse(((eye(size(Pb_01))-Pb_01)^-1));

v_0 = zeros(1,size(Pb_01,1));
v_0(1:size(Pb_10,1)) = xbc*Pbc_b_10/(xbc*Pbc_b_10*ones(size(Pbc_b_10,2),1));

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


%%%%%%%%%%%%%%%%%% centered moments residence time in state 1 %%%%%%%%%%%%%%%%%%%%%%%%%
xb = x(1:switch_indx_01-1);
xbc = x(switch_indx_01:end);

v_1 = zeros(1,size(Pbc_10,1));
v_1(switch_indx_01-size(Pb_10,1):end) = xb*Pb_bc_01/(xb*Pb_bc_01*ones(size(Pb_bc_01,2),1));

% v_1 = xb*Pb_bc_10/(xb*Pb_bc_10*ones(size(Pb_bc_10,2),1));

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