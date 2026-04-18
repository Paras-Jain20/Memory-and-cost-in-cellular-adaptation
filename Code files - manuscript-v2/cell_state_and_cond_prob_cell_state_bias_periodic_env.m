function [prob_state_1,prob_state_1_given_1,mean_pi,P_Sn_En,P_s] = cell_state_and_cond_prob_cell_state_bias_periodic_env(N,T,q01,q10,Eh_bias,varargin)

if(~isempty(varargin))
    fidelity_00 = varargin{1}(1);
    fidelity_11 = varargin{1}(2);
else
    fidelity_00 = 1;
    fidelity_11 = 1;
end

El_bias = 1-Eh_bias;
T_El = floor(T*(1-Eh_bias));
T_Eh = T - T_El;

P0 = transition_matrix_P_cal(N,1-fidelity_00);
P1 = transition_matrix_P_cal(N,fidelity_11);

%% stationary distribution calculation
prob_analy = zeros(T,N+1);
prob_analy(end,:) = null((eye(N+1) - P0^(T_El)*P1^(T_Eh))');
prob_analy(end,:)  = prob_analy(end,:)/sum(prob_analy(end,:));
mean_pi = zeros(1,T);
for indx = 1:T-1
    if(indx <= T_El)
        prob_analy(indx,:) = prob_analy(end,:)*P0^(indx);
    else
        prob_analy(indx,:) = prob_analy(end,:)*P0^(T_El)*P1^(indx-T_El);
    end
    mean_pi(indx) = sum((0:N).*prob_analy(indx,:))/N;
end
mean_pi(end) = sum((0:N).*prob_analy(end,:))/N;

%% Partitioning the state transition matrices
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

x_0 = zeros(total_time_steps+1,size(P0b_01,1));
x_1 = zeros(total_time_steps,size(P0bc_10,1));

if(sum(sum(xbc))~=0)
    x_1(1,:) = xbc/sum(xbc);
end
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
    
        if(rem(k-l*T,T) < T_El)
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
            if(sum(abs(x_0(k+1,:)-x_0(k+1-T,:)))< 10^-4 && sum(abs(x_1(k+1,:)-x_1(k+1-T,:)))< 10^-4 && rem(k+1,T) == 0)
                temp_prob_state_1(1,k+1) = sum(x_1(k+1,:));
                break;
            end
        end
        k = k+1;
    end

% temp_prob_state_1 = temp_prob_state_1(temp_prob_state_1=0,:);
% x_1 = x_1(~all(x_1==0,2),:);
x_1 = x_1(1:k+1,:);
x_0 = x_0(1:k+1,:);
prob_state_1 = temp_prob_state_1(k+1-T+1:k+1);

%% Calculating P(Sm = 1, Sm-1 = 1)
prob_state_1_given_1 = zeros(1,T);
for m = 1:T
    temp_x_0 = zeros(size(x_0(1,:)));
    temp_x_1 = x_1(end-T+m-1,:);

    if(sum(temp_x_1) ==0)
        prob_state_1_given_1(m) = sum(temp_x_1);
    else
        temp_x_1 = temp_x_1/sum(temp_x_1);
        if(m <= T_El)
            temp = temp_x_1*P0bc_10;
            temp_x_1 = [temp(1:floor(q01*N)+1-size(P0b_10,1)) temp_x_0*P0b_bc_01+temp(floor(q01*N)+2-size(P0b_10,1):end)];
        else
            temp = temp_x_1*P1bc_10;
            temp_x_1 = [temp(1:floor(q01*N)+1-size(P1b_10,1)) temp_x_0*P1b_bc_01+temp(floor(q01*N)+2-size(P1b_10,1):end)];
        end
        prob_state_1_given_1(m) = sum(temp_x_1);
    end
end

%% Calculation of P_Sn_En 
temp_prob_state_1 = prob_state_1([T 1:T-1]);

P_Sn_En = zeros(T,4);
env = [zeros(1,T_El) ones(1,T_Eh)];
for m = 1:T
    if(m == 1)
        P_Sn_En(m,1) = env(T) * env(m) * temp_prob_state_1(m) + (1-env(T)) * env(m) * temp_prob_state_1(m); % P_Eh_Sh
        P_Sn_En(m,2) = env(T) * (1-env(m)) * temp_prob_state_1(m) + (1-env(T)) * (1-env(m)) * temp_prob_state_1(m); % P_El_Sh
        P_Sn_En(m,3) = env(T) * env(m) * (1-temp_prob_state_1(m)) + (1-env(T)) * env(m) * (1-temp_prob_state_1(m)); % P_Eh_Sl
        P_Sn_En(m,4) = env(T) * (1-env(m)) * (1-temp_prob_state_1(m)) + (1-env(T)) * (1-env(m)) * (1-temp_prob_state_1(m)); % P_El_Sl
    else
        P_Sn_En(m,1) = env(m-1) * env(m) * temp_prob_state_1(m) + (1-env(m-1)) * env(m) * temp_prob_state_1(m); % P_Eh_Sh
        P_Sn_En(m,2) = env(m-1) * (1-env(m)) * temp_prob_state_1(m) + (1-env(m-1)) * (1-env(m)) * temp_prob_state_1(m); % P_El_Sh
        P_Sn_En(m,3) = env(m-1) * env(m) * (1-temp_prob_state_1(m)) + (1-env(m-1)) * env(m) * (1-temp_prob_state_1(m)); % P_Eh_Sl
        P_Sn_En(m,4) = env(m-1) * (1-env(m)) * (1-temp_prob_state_1(m)) + (1-env(m-1)) * (1-env(m)) * (1-temp_prob_state_1(m)); % P_El_Sl
    end
end

%% Determining transition probabilities between phenotypic states

%P(S_n+1 = S_h|S_n = S_l) = P(k_n+1 = floor(q01*N)+1 |k_n = floor(q01*N), S_n = S_l) P(k_n = floor(q01*N)|S_n = S_l)
%P(S_n+1 = S_l|S_n = S_h) = P(k_n+1 = floor(q10*N) |k_n = floor(q10*N)+1, S_n = S_h) P(k_n = floor(q10*N)+1|S_n = S_h)
P_s = zeros(2,2,T);
for m = 1:T
    temp_x_0 = x_0(end-T+m-2,:);
    temp_x_1 = x_1(end-T+m-2,:);

    if(m==1)
        if(sum(temp_x_0) == 0)
            t_0_1 = 0;
            t_0_0 = 1-t_0_1;
        else
            t_0_1 = P1((floor(q01*N)+1),(floor(q01*N)+1)+1)*temp_x_0(end)/sum(temp_x_0);
            t_0_0 = 1-t_0_1;
        end

        if(sum(temp_x_1) == 0)
            t_1_0 = 0;
            t_1_1 = 1-t_1_0;
        else
            t_1_0 = P1((floor(q10*N)+1)+1,(floor(q10*N)+1))*temp_x_1(1)/sum(temp_x_1);
            t_1_1 = 1-t_1_0;
        end
    elseif(m <= T_El+1)
        if(sum(temp_x_0) == 0)
            t_0_1 = 0;
            t_0_0 = 1-t_0_1;
        else
            t_0_1 = P0((floor(q01*N)+1),(floor(q01*N)+1)+1)*temp_x_0(end)/sum(temp_x_0);
            t_0_0 = 1-t_0_1;
        end

        if(sum(temp_x_1) == 0)
            t_1_0 = 0;
            t_1_1 = 1-t_1_0;
        else
            t_1_0 = P0((floor(q10*N)+1)+1,(floor(q10*N)+1))*temp_x_1(1)/sum(temp_x_1);
            t_1_1 = 1-t_1_0;
        end
    else
        if(sum(temp_x_0) == 0)
            t_0_1 = 0;
            t_0_0 = 1-t_0_1;
        else
            t_0_1 = P1((floor(q01*N)+1),(floor(q01*N)+1)+1)*temp_x_0(end)/sum(temp_x_0);
            t_0_0 = 1-t_0_1;
        end

        if(sum(temp_x_1) == 0)
            t_1_0 = 0;
            t_1_1 = 1-t_1_0;
        else
            t_1_0 = P1((floor(q10*N)+1)+1,(floor(q10*N)+1))*temp_x_1(1)/sum(temp_x_1);
            t_1_1 = 1-t_1_0;
        end
    end
    P_s(:,:,m) = [t_0_0 t_0_1; t_1_0 t_1_1];
end

end