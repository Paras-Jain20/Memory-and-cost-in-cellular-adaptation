function centered_moments_adapt_time = undated_memory_adaptation_time_periodic_env(N,T,q01,q10)
centered_moments_adapt_time = zeros(2,2);

P0 = transition_matrix_P_cal(N,0);
P1 = transition_matrix_P_cal(N,1);

% stationary distribution calculation
prob_analy = zeros(T,N+1);
prob_analy(end,:) = null((eye(N+1) - P0^(floor(T/2))*P1^(ceil(T/2)))');
prob_analy(end,:)  = prob_analy(end,:)/sum(prob_analy(end,:));
for indx = 1:T-1
    if(indx <= floor(T/2))
        prob_analy(indx,:) = prob_analy(end,:)*P0^(indx);
    else
        prob_analy(indx,:) = prob_analy(end,:)*P0^(floor(T/2))*P1^(indx-floor(T/2));
    end
end

%%% Partitioning the state transition matrices
P0b_01 = P0(1:floor(q01*N)+1,1:floor(q01*N)+1);
P0b_bc_01 = P0(1:floor(q01*N)+1,floor(q01*N)+2:end);
P0bc_b_01 = P0(floor(q01*N)+2:end,1:floor(q01*N)+1);
P0bc_01 = P0(floor(q01*N)+2:end,floor(q01*N)+2:end);

P1b_01 = P1(1:floor(q01*N)+1,1:floor(q01*N)+1);
P1b_bc_01 = P1(1:floor(q01*N)+1,floor(q01*N)+2:end);
P1bc_b_01 = P1(floor(q01*N)+2:end,1:floor(q01*N)+1);
P1bc_01 = P1(floor(q01*N)+2:end,floor(q01*N)+2:end);

P0b_10 = P0(1:floor(q10*N)+1,1:floor(q10*N)+1);
P0b_bc_10 = P0(1:floor(q10*N)+1,floor(q10*N)+2:end);
P0bc_b_10 = P0(floor(q10*N)+2:end,1:floor(q10*N)+1);
P0bc_10 = P0(floor(q10*N)+2:end,floor(q10*N)+2:end);

P1b_10 = P1(1:floor(q10*N)+1,1:floor(q10*N)+1);
P1b_bc_10 = P1(1:floor(q10*N)+1,floor(q10*N)+2:end);
P1bc_b_10 = P1(floor(q10*N)+2:end,1:floor(q10*N)+1);
P1bc_10 = P1(floor(q10*N)+2:end,floor(q10*N)+2:end);

%% Calculating P(Sm = 1), where m = [1,T]

xbc = prob_analy(end,floor(q10*N)+2:end);

n = 10; % run the simulation for n*T timestep to reach saturation
total_time_steps = n*T;

x_0 = zeros(total_time_steps,size(P0b_01,1));
x_1 = zeros(total_time_steps,size(P0bc_10,1));

x_1(1,:) = xbc/sum(xbc);
% x_1(1,floor(q01*N)+2-size(Pb_10,1):end) = xbc;

temp_prob_state_1 = zeros(1, total_time_steps);
temp_prob_state_0 = zeros(1, total_time_steps);
l = 0;
k = 1;
    while(true)
    
        temp_prob_state_1(1,k) = sum(x_1(k,:));
        temp_prob_state_0(1,k) = sum(x_0(k,:));
    
        if(mod(k,T) == 0)
            l = l+1;
        end
    
        if(rem(k-l*T,T) < floor(T/2))
            temp = x_1(k,:)*P0bc_10;
            x_1(k+1,:) = [temp(1:floor(q01*N)+1-size(P0b_10,1)) x_0(k,:)*P0b_bc_01+temp(floor(q01*N)+2-size(P0b_10,1):end)];
    
            temp = x_0(k,:)*P0b_01;
            x_0(k+1,:) = [temp(1:floor(N*q10)+1)+x_1(k,:)*P0bc_b_10 temp(floor(N*q10)+2:end)];
        else
            temp = x_1(k,:)*P1bc_10;
            x_1(k+1,:) = [temp(1:floor(q01*N)+1-size(P1b_10,1)) x_0(k,:)*P1b_bc_01+temp(floor(q01*N)+2-size(P1b_10,1):end)];
    
            temp = x_0(k,:)*P1b_01;
            x_0(k+1,:) = [temp(1:floor(N*q10)+1)+x_1(k,:)*P1bc_b_10 temp(floor(N*q10)+2:end)];
        end
        if(k+1> T)
            if(sum(abs(x_0(k+1,:)-x_0(k+1-T,:)))< 10^-6 && rem(k+1,T) == 0)
                temp_prob_state_1(1,k+1) = sum(x_1(k+1,:));
                break;
            end
        end
        k = k+1;
    end

% temp_prob_state_1 = temp_prob_state_1(temp_prob_state_1=0,:);
x_1 = x_1(~all(x_1==0,2),:);
x_0 = x_0(~all(x_1==0,2),:);

x_1_ss = x_1(end-T+1:end,:);
x_0_ss = x_0(end-T+1:end,:);

%% calculating time to adaptation in a periodic environment : adaptation from state 0 to 1 when sequence of 1's are encountered

v_0 = x_0_ss(floor(T/2),:)/sum(x_0_ss(floor(T/2),:));
% v_0_temp = prob_analy(floor(T/2),1:floor(q01*N)+1)/(sum(prob_analy(floor(T/2),1:floor(q01*N)+1)));

total_time_steps = 10000;
prob_adapt_time_01 = zeros(1, total_time_steps);
mean_adapt_time_01 = 0;
expect_adapt_time_sqr_01 = 0;
m = 0;
for k = 1:total_time_steps
    if(mod(k,T) == 0)
        m = m+1;
    end
    if(k< ceil(T/2))
        prob_adapt_time_01(k) = sum(v_0 * P1b_01^(k-1) * P1b_bc_01);
    elseif(k >= ceil(T/2) && k < T)
        prob_adapt_time_01(k) = sum(v_0 * P1b_01^(ceil(T/2)-1) * P0b_01^(k-ceil(T/2)) * P0b_bc_01);
    elseif(k >= m*T && k < m*T + ceil(T/2) )
        prob_adapt_time_01(k) = sum(v_0 * P1b_01^(ceil(T/2)-1) * P0b_01^(floor(T/2)) *(P1b_01^(ceil(T/2)) * P0b_01^(floor(T/2)))^(m-1) * P1b_01^(k-m*T) * P1b_bc_01);
    elseif(k >= m*T + ceil(T/2)  && k < (m+1)*T)
        prob_adapt_time_01(k) = sum(v_0 * P1b_01^(ceil(T/2)-1) * P0b_01^(floor(T/2)) *(P1b_01^(ceil(T/2)) * P0b_01^(floor(T/2)))^(m-1) * P1b_01^(ceil(T/2)) * P0b_01^(k-(m*T+ceil(T/2))) * P0b_bc_01);
    end
    mean_adapt_time_01 = mean_adapt_time_01 + k*prob_adapt_time_01(k);
    expect_adapt_time_sqr_01 = expect_adapt_time_sqr_01 + k^2 *prob_adapt_time_01(k);
end
var_adapt_time_01 = expect_adapt_time_sqr_01 - mean_adapt_time_01^2;

centered_moments_adapt_time(:,1) = [mean_adapt_time_01;var_adapt_time_01];

% figure
% bar(prob_adapt_time_01(1:100));

%% calculating time to adaptation in a periodic environment : adaptation from state 1 to 0 when sequence of 0's are encountered
v_1 = x_1_ss(end,:)/sum(x_1_ss(end,:));
% v_1_temp = prob_analy(end,floor(q10*N)+2:end)/(sum(prob_analy(end,floor(q10*N)+2:end)));

total_time_steps = 10000;

% another way of plotting residence time distribution
prob_adapt_time_10 = zeros(1, total_time_steps);
mean_adapt_time_10 = 0;
expect_adapt_time_sqr_10 = 0;

m = 0;
for k = 1:total_time_steps
    if(mod(k,T) == 0)
        m = m+1;
    end
    if(k< floor(T/2))
        prob_adapt_time_10(k) = sum(v_1 * P0bc_10^(k-1) * P0bc_b_10);
    elseif(k >= floor(T/2) && k < T)
        prob_adapt_time_10(k) = sum(v_1 * P0bc_10^(floor(T/2)-1) * P1bc_10^(k-floor(T/2)) * P1bc_b_10);
    elseif(k >= m*T && k < m*T + floor(T/2) )
        prob_adapt_time_10(k) = sum(v_1 * P0bc_10^(floor(T/2)-1) * P1bc_10^(ceil(T/2)) *(P0bc_10^(floor(T/2)) * P1bc_10^(ceil(T/2)))^(m-1) * P0bc_10^(k-m*T) * P0bc_b_10);
    elseif(k >= m*T + floor(T/2)  && k < (m+1)*T)
        prob_adapt_time_10(k) = sum(v_1 * P0bc_10^(floor(T/2)-1) * P1bc_10^(ceil(T/2)) *(P0bc_10^(floor(T/2)) * P1bc_10^(ceil(T/2)))^(m-1) * P0bc_10^(floor(T/2)) * P1bc_10^(k-(m*T+floor(T/2))) * P1bc_b_10);
    end
    mean_adapt_time_10 = mean_adapt_time_10 + k*prob_adapt_time_10(k);
    expect_adapt_time_sqr_10 = expect_adapt_time_sqr_10 + k^2 *prob_adapt_time_10(k);
end
var_adapt_time_10 = expect_adapt_time_sqr_10 - mean_adapt_time_10^2;

centered_moments_adapt_time(:,2) = [mean_adapt_time_10;var_adapt_time_10];
end