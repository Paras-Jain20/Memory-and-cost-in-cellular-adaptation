function [centered_moments] = undated_memory_marginal_resi_time_moments(N,q01,q10,P)
centered_moments = zeros(4,2);
uncentered_moments = zeros(4,2);

%% Stationary sitribution calculations
p = P(1,2);
x = zeros(1,N+1);
for x_indx = 1:N+1
    x(x_indx) = nchoosek(N,x_indx-1)*(p)^(x_indx-1)* (1-p)^(N-(x_indx-1));
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

xb = x(1:floor(q10*N)+1);
xbc = x(floor(q10*N)+2:end);
I_Pb_01_inv = sparse(((eye(size(Pb_01))-Pb_01)^-1));

%% Residence time moment calculation

%%%%%%%%%%%%%%%%%% centered moments residence time in state 0 %%%%%%%%%%%%%%%%%%%%%%%%%

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

xb = x(1:floor(q01*N)+1);
xbc = x(floor(q01*N)+2:end);

v_1 = zeros(1,size(Pbc_10,1));
v_1(floor(q01*N)+2-size(Pb_10,1):end) = xb*Pb_bc_01/(xb*Pb_bc_01*ones(size(Pb_bc_01,2),1));

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