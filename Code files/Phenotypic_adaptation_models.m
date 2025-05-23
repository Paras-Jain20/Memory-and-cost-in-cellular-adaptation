%%%% Simulation dynamics of cell state with memory models
function out = Phenotypic_adaptation_models
out{1} = @Costly_transition;
out{2}= @Markovian_Costly_transition;

function [cell_state,benefit,avg_benefit,avg_benefit_with_time,pi_est] = Costly_transition(reward,cost,total_time_steps, env, N, p, c0)
% env_memory_store = zeros(total_time_steps+1,N);
cell_state = zeros(total_time_steps,1); % cell states corresponding to each environment and time point
env_memory = flip(env(1:N)); %binornd(1,p,[N,1]); % assumming that the past N environment encounter are all 0 state
% env_memory_store(1,:) = env_memory;

benefit = zeros(total_time_steps,1);
env_inf = mean(env_memory); % inferred environmental state at time step 1
pi_est = zeros(total_time_steps,1);
avg_benefit_with_time = zeros(total_time_steps,1);

% to calculate the value at the starting timepoint
q = (reward(2) - cost(1))/(sum(reward) - sum(cost));% indifference probability when there is no cost associated with transition
q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

if(env_inf <= q)
    cell_state(1) = 0;
elseif(env_inf > q)
    cell_state(1) = 1;
end

env = env(N+1:end);

if(cell_state(1) == 0 && env(1) == 0)
    benefit(1) = benefit(1) + reward(2);
elseif(cell_state(1) == 0 && env(1) == 1)
    benefit(1) = benefit(1) + cost(2);
elseif(cell_state(1) == 1 && env(1) == 1)
    benefit(1) = benefit(1) + reward(1);
elseif(cell_state(1) == 1 && env(1) == 0)
    benefit(1) = benefit(1) + cost(1);
end

env_memory(2:end,:) = env_memory(1:end-1,:);
env_memory(1,:) = env(1,:);
env_inf = mean(env_memory);

for t = 1:total_time_steps-1 % t keep track of time steps
        % q = (reward(env_axis,2) - cost(env_axis,1))/(sum(reward) - sum(cost));% indifference probability when there is no cost associated with transition

        % below we are not anticipating the next time step environment

        if(cell_state(t)  == 0 && env_inf <= q01)
            cell_state(t+1) = 0;
        elseif(cell_state(t)  == 0 && env_inf > q01)
            cell_state(t+1) = 1;
            benefit(t+1) = -c0;
        elseif(cell_state(t)  == 1 && env_inf > q10)
            cell_state(t+1) = 1; % paying the cost of transition
        elseif(cell_state(t)  == 1 && env_inf <= q10)
            cell_state(t+1) = 0;
            benefit(t+1) = -c0; % paying the cost of transition
        end

        if(cell_state(t+1) == 0 && env(t+1) == 0)
            benefit(t+1) = benefit(t+1) + reward(2);
        elseif(cell_state(t+1) == 0 && env(t+1) == 1)
            benefit(t+1) = benefit(t+1) + cost(2);
        elseif(cell_state(t+1) == 1 && env(t+1) == 1)
            benefit(t+1) = benefit(t+1) + reward(1);
        elseif(cell_state(t+1) == 1 && env(t+1) == 0)
            benefit(t+1) = benefit(t+1) + cost(1);
        end

        avg_benefit_with_time(t,:) = reward(1) * env_inf + cost(1) * (1-env_inf);
        
%         if(cell_state(t+1,env_axis) == 1)
%         avg_benefit_with_time(t,:) = reward(env_axis,1) * env_inf(env_axis) + cost(env_axis,1) * (1-env_inf(env_axis));
%         elseif(cell_state(t+1,env_axis)  == 0)
%         avg_benefit_with_time(t,:) = reward(env_axis,2) * (1-env_inf(env_axis)) + cost(env_axis,2) * (env_inf(env_axis));
%         end

    
    % env_memory_store(t+1,:) = env_memory;
    pi_est(t,:) = env_inf;
    env_memory(2:end,:) = env_memory(1:end-1,:);
    env_memory(1,:) = env(t+1,:);
    env_inf = mean(env_memory);
     
end
avg_benefit = mean(benefit(100:end,:),1);

function [cell_state,benefit,avg_benefit,avg_benefit_with_time,pi_est] = Markovian_Costly_transition(reward,cost,total_time_steps, env, N, p, c0)

cell_state = zeros(total_time_steps,1); % cell states corresponding to each environment and time point

benefit = zeros(total_time_steps,1);
k_inf = sum(env(1:N)); % inferred environmental state at time step 1
pi_est = zeros(total_time_steps,1);
avg_benefit_with_time = zeros(total_time_steps,1);

% to calculate the value at the starting timepoint
% q = (reward(2) - cost(1))/(sum(reward) - sum(cost));% indifference probability when there is no cost associated with transition
q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

if(cell_state(1) == 0 && k_inf <= q01*N)
    cell_state(1) = 0;
elseif(cell_state(1) == 0 && k_inf > q01*N)
    cell_state(1) = 1;
    benefit(1) = -c0;
elseif(cell_state(1) == 1 && k_inf > q10*N)
    cell_state(1) = 1;    
elseif(cell_state(1) == 1 && k_inf <= q10*N)
    cell_state(1) = 0;
    benefit(1) = -c0;
end

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

reaction_prob = [1-k_inf/N k_inf/N];

for i = 1:length(reaction_prob)
    if (sum(reaction_prob(1:i)) >= rand_set(1))
        react_type = i;
        break;
    end
end

if(react_type == 1)
    k_inf = k_inf + env(1);
else
    k_inf = k_inf - 1 + env(1);
end

for t = 1:total_time_steps-1 % t keep track of time steps
    % q = (reward(env_axis,2) - cost(env_axis,1))/(sum(reward) - sum(cost));% indifference probability when there is no cost associated with transition

    % below we are not anticipating the next time step environment


    if(cell_state(t) == 0 && k_inf <= q01*N)
        cell_state(t+1) = 0;
    elseif(cell_state(t) == 0 && k_inf > q01*N)
        cell_state(t+1) = 1;
        benefit(t+1) = -c0;
    elseif(cell_state(t) == 1 && k_inf > q10*N)
        cell_state(t+1) = 1;
    elseif(cell_state(t) == 1 && k_inf <= q10*N)
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

    if(react_type == 1)
        k_inf = k_inf + env(t+1);
    else
        k_inf = k_inf - 1 + env(t+1);
    end

end
avg_benefit = mean(benefit(100:end,:),1);