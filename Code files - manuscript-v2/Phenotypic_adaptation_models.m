function [cell_state,benefit,avg_benefit,avg_benefit_with_time,pi_est] = Phenotypic_adaptation_models(reward,cost,total_time_steps, env, N, p, c0,varargin)

if(~isempty(varargin))
    fidelity_00 = varargin{1}(1);
    fidelity_11 = varargin{1}(2);
else
    fidelity_00 = 1;
    fidelity_11 = 1;
end

cell_state = zeros(total_time_steps,1); % cell states corresponding to each environment and time point

benefit = zeros(total_time_steps,1);
k_inf = sum(env(1:N)); % inferred environmental state at time step 1
pi_est = zeros(total_time_steps,1);
avg_benefit_with_time = zeros(total_time_steps,1);

% to calculate the value at the starting timepoint
% q = (reward(2) - cost(1))/(sum(reward) - sum(cost));% indifference probability when there is no cost associated with transition
q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

% if(cell_state(1) == 0 && k_inf <= q01*N)
%     cell_state(1) = 0;
% elseif(cell_state(1) == 0 && k_inf > q01*N)
%     cell_state(1) = 1;
%     benefit(1) = -c0;
% elseif(cell_state(1) == 1 && k_inf > floor(q10*N))
%     cell_state(1) = 1;    
% elseif(cell_state(1) == 1 && k_inf <= floor(q10*N))
%     cell_state(1) = 0;
%     benefit(1) = -c0;
% end

env = env(N+1:end);

if(cell_state(1) == 0 && env(1) == 0)
    benefit(1) = benefit(1)+ reward(2);
elseif(cell_state(1) == 0 && env(1) == 1)
    benefit(1) = benefit(1)+cost(2);
elseif(cell_state(1) == 1 && env(1) == 1)
    benefit(1) = benefit(1)+reward(1);
elseif(cell_state(1) == 1 && env(1) == 0)
    benefit(1) = benefit(1)+cost(1);
end

rand_set = rand(total_time_steps+1,1);
rand_set_1 = rand(total_time_steps+1,1);

reaction_prob = [1-k_inf/N k_inf/N];

for i = 1:length(reaction_prob)
    if (sum(reaction_prob(1:i)) >= rand_set(1))
        react_type = i;
        break;
    end
end

if(env(1) == 0)
    if(fidelity_00 > rand_set_1(1))
        temp = 0;
    else
        temp = 1;
    end
else
    if(fidelity_11 > rand_set_1(1))
        temp = 1;
    else
        temp = 0;
    end
end
if(react_type == 1)
    k_inf = k_inf + temp;
else
    k_inf = k_inf - 1 + temp;
end

for t = 1:total_time_steps-1 % t keep track of time steps
    % q = (reward(env_axis,2) - cost(env_axis,1))/(sum(reward) - sum(cost));% indifference probability when there is no cost associated with transition

    % below we are not anticipating the next time step environment

    if(cell_state(t) == 0 && k_inf <= floor(q01*N))
        cell_state(t+1) = 0;
    elseif(cell_state(t) == 0 && k_inf > floor(q01*N))
        cell_state(t+1) = 1;
        benefit(t+1) = -c0;
    elseif(cell_state(t) == 1 && k_inf > floor(q10*N))
        cell_state(t+1) = 1;
    elseif(cell_state(t) == 1 && k_inf <= floor(q10*N))
        cell_state(t+1) = 0;
        benefit(t+1) = -c0;
    end


    if(cell_state(t+1) == 0 && env(t+1) == 0)
        benefit(t+1) = benefit(t+1)+reward(2);
    elseif(cell_state(t+1) == 0 && env(t+1) == 1)
        benefit(t+1) = benefit(t+1)+cost(2);
    elseif(cell_state(t+1) == 1 && env(t+1) == 1)
        benefit(t+1) = benefit(t+1)+reward(1);
    elseif(cell_state(t+1) == 1 && env(t+1) == 0)
        benefit(t+1) = benefit(t+1)+cost(1);
    end

    avg_benefit_with_time(t,:) = reward(1) * k_inf/N + cost(1) * (1-k_inf/N);

    pi_est(t,:) = k_inf/N;

    reaction_prob = [1-k_inf/N k_inf/N];

    for i = 1:length(reaction_prob)
        if (sum(reaction_prob(1:i)) >= rand_set(t+1))
            react_type = i;
            break;
        end
    end

    if(env(t+1) == 0)
        if(fidelity_00 > rand_set_1(t+1))
            temp = 0;
        else
            temp = 1;
        end
    else
        if(fidelity_11 > rand_set_1(t+1))
            temp = 1;
        else
            temp = 0;
        end
    end
    if(react_type == 1)
        k_inf = k_inf + temp;
    else
        k_inf = k_inf - 1 + temp;
    end
end
avg_benefit = mean(benefit(100:end,:),1);
end