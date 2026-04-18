analysis_time = [24 24*15];
total_time_steps = sum(analysis_time);

N = 40;%[4 8 12]; %11; % memory capacity

cost = [0 -1].*(1/20);
reward = [0.2 1].*(1/20);
% cost = [0.3 0.1].*(1/20);
% reward = [0.8 1].*(1/20);

q = (reward(2)-cost(1))/(sum(reward) - sum(cost));
alpha = reward(1)-cost(2);
beta = alpha*q./(1-q); % indifference environment (p_I)

one_step_auto_corr = 0;

mean_env_pair = [0 0.2;
    0 0.4;
    0 0.6;
    0 0.8;
    0 1];

% mean_env_pair = [zeros(60,1) [0.4:0.01:0.99]'];

c0_set = 0;%linspace(0,min(alpha,beta)-min(alpha,beta)/100,3);

%% Analytical Calculations: One-step shift in the environment and single clonal population 
adaptation_time = zeros(length(c0_set),size(mean_env_pair,1)); % first two indices corresponds to the number of p and model values
benefit_time_series = zeros(analysis_time(2),length(c0_set),size(mean_env_pair,1)); % first two indices corresponds to the number of p and model values
cell_state_time_series = zeros(analysis_time(2),length(c0_set),size(mean_env_pair,1)); % first two indices corresponds to the number of p and model values

c0_length = length(c0_set);
mean_env_pair_length = size(mean_env_pair,1);

for c0_indx = 1:c0_length
    q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
    q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
    for mean_env_pair_indx = 1:mean_env_pair_length
        disp(['running analysis for c0 ' num2str(c0_indx) ' mean_env_pair ' num2str(mean_env_pair_indx) ])
        p = mean_env_pair(mean_env_pair_indx,2);

        [centered_moments,~,prob_state_1_start_0,prob_Sn_Sh_Sn_1_Sh_start_0,prob_Sn_Sl_Sn_1_Sl_start_0] = one_step_env_switch_mean_adaptation_time(N,q01,q10,mean_env_pair(mean_env_pair_indx,1),mean_env_pair(mean_env_pair_indx,2),analysis_time(2));
        adaptation_time(c0_indx,mean_env_pair_indx) = centered_moments(1,1);
        
        P_s = zeros(2,2,length(prob_Sn_Sl_Sn_1_Sl_start_0));
        P_s(1,1,:) = prob_Sn_Sl_Sn_1_Sl_start_0;
        P_s(1,2,:) = 1-P_s(1,1,:);
        P_s(2,2,:) = prob_Sn_Sh_Sn_1_Sh_start_0;
        P_s(2,1,:) = 1-P_s(2,2,:);
        
        % Dynamics of population cell fraction and cellular benefit
        x = zeros(analysis_time(2),3);
        fd = zeros(analysis_time(2),1);
        x([1 2],:) = repmat([1 0 0],2,1);
        benefit_time_series(1,c0_indx,mean_env_pair_indx) = reward(2);

        for n = 2:analysis_time(2)-1
            if((1-prob_state_1_start_0(n))~=0)
                S_l_fit = reward(2)*(1-p) +cost(2)*(p)-prob_state_1_start_0(n-1)*P_s(2,1,n)*c0_set(c0_indx)/(1-prob_state_1_start_0(n));
            else
                S_l_fit = reward(2)*(1-p) +cost(2)*(p);
            end

            if((prob_state_1_start_0(n))~=0)
                S_h_fit = reward(1)*(p)+cost(1)*(1-p)-(1-prob_state_1_start_0(n-1))*P_s(1,2,n)*c0_set(c0_indx)/prob_state_1_start_0(n);
            else
                S_h_fit = reward(1)*(p)+cost(1)*(1-p);
            end

            transiton_matrix = zeros(2+1);
            growth_matrix = zeros(2+1);

            transiton_matrix(1:2, 3) = [-S_l_fit*(1-double(S_l_fit>0)); -S_h_fit*(1-double(S_h_fit>0))];
            transiton_matrix(1:2,1:2) = (1-transiton_matrix(1:2,3)).*[P_s(1,1,n) P_s(1,2,n); P_s(2,1,n) P_s(2,2,n)];
            growth_matrix(1:2,1:2) = [S_l_fit*double(S_l_fit>0) 0; 0 S_h_fit*double(S_h_fit>0)];

            transiton_matrix(3,:) = [0 0 1];

            x(n+1,:) = x(n,:)*(transiton_matrix+growth_matrix);
            fd(n) = x(n,:)*sum(transiton_matrix+growth_matrix,2);
            x(n+1,:) = x(n+1,:)/fd(n);

            benefit_time_series(n,c0_indx,mean_env_pair_indx) = S_l_fit*(x(n,1)./sum(x(n,[1 2]),2)) + S_h_fit*(x(n,2)./sum(x(n,[1 2]),2));
        end
            cell_state_time_series(:,c0_indx,mean_env_pair_indx) = x(:,2)./sum(x(:,[1 2]),2);
    end
end
%% Numerical Calculations: One-step shift in the environment and single clonal population 
num_rep= 10000; % number of replicate for numberical simulations

num_benefit_centered_moments = zeros(4,length(c0_set),size(mean_env_pair,1),num_rep); % first two indices corresponds to the number of p and model values
num_benefit_time_series = zeros(total_time_steps,length(c0_set),size(mean_env_pair,1),num_rep); % first two indices corresponds to the number of p and model values
num_cell_state_time_series = zeros(total_time_steps,length(c0_set),size(mean_env_pair,1),num_rep); % first two indices corresponds to the number of p and model values
env_time_series = zeros(total_time_steps,length(c0_set),size(mean_env_pair,1),num_rep);
c0_length = length(c0_set);
mean_env_pair_length = size(mean_env_pair,1);

for c0_indx = 1:c0_length
    q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
    q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

    parfor mean_env_pair_indx = 1:mean_env_pair_length
        for rep_indx = 1:num_rep
            disp(['running analysis for rep num ' num2str(rep_indx) ' c0 ' num2str(c0_indx) ' mean_env_pair ' num2str(mean_env_pair_indx) ])

            p = mean_env_pair(mean_env_pair_indx,1);
            p2 = min(1,(1-one_step_auto_corr)*(1-p));
            p1 = min(1,(1-one_step_auto_corr)*(p));

            env_1 = correlated_env(p1,p2,analysis_time(1)+N); % the state of the environment along each axis and timesteps

            p = mean_env_pair(mean_env_pair_indx,2);
            p2 = min(1,(1-one_step_auto_corr)*(1-p));
            p1 = min(1,(1-one_step_auto_corr)*(p));

            env_2 = correlated_env(p1,p2,analysis_time(2)); % the state of the environment along each axis and timesteps

            env = [env_1 env_2]';

            [cell_state_with_time,benefit,avg_benefit,avg_benefit_with_time,~] = Phenotypic_adaptation_models(reward,cost,total_time_steps, env, N, p, c0_set(c0_indx));
            temp_auto_corr = autocorr(benefit);

            num_benefit_centered_moments(:,c0_indx,mean_env_pair_indx,rep_indx) = [avg_benefit; var(benefit); skewness(benefit); temp_auto_corr(2)]; % saving lag 1 autocorrelation
            num_benefit_time_series(:,c0_indx,mean_env_pair_indx,rep_indx) = benefit;
            num_cell_state_time_series(:,c0_indx,mean_env_pair_indx,rep_indx) = cell_state_with_time;
            env_time_series(:,c0_indx,mean_env_pair_indx,rep_indx) = env(N+1:end);
            temp_adaptation_time = find(cell_state_with_time == 1,1)-total_time_steps;
            if(~isempty(temp_adaptation_time))
                adaptation_time_num(c0_indx,mean_env_pair_indx,rep_indx) = find(cell_state_with_time==1,1)-total_time_steps-1;
            else
                adaptation_time_num(c0_indx,mean_env_pair_indx,rep_indx) = nan;
            end
        end
    end
end
%% single clonal population 
%%%%%%%%%%%%%%% Cell state fraction and per-capita growth dynamics with time; variable mean_env_pair, fixed cost %%%%%%%%%%%

for out_type = 1:2 % 1: cell state, 2: avg benefit
    figure('Position',[597 568 430 236])
    color_order = ["#7E2F8E", "#D95319","#EDB120","#77AC30","#0072BD"];
    upper_value =  24*10;
    lower_value = 24;
    for mean_env_pair_indx = size(mean_env_pair,1):-1:1
        for c0_indx = 1%:2
            if(out_type == 1)
                % data = reshape(num_cell_state_time_series(:,c0_indx,mean_env_pair_indx,:),2*total_time_steps,num_rep);
                y_analy = [zeros(analysis_time(1),1); cell_state_time_series(2:upper_value+2,c0_indx,mean_env_pair_indx)];
            else
                % data = reshape(num_benefit_time_series(:,c0_indx,mean_env_pair_indx,:),2*total_time_steps,num_rep);
                y_analy = [ones(analysis_time(1),1)*reward(2); benefit_time_series(2:upper_value+2,c0_indx,mean_env_pair_indx)];

            end
            % y = mean(data(total_time_steps-lower_value+1:total_time_steps+upper_value+1,:),2);
            plot(analysis_time(1)-lower_value:analysis_time(1)+upper_value,y_analy(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1),'LineWidth',2,Color=color_order(mean_env_pair_indx));
            hold on
            % plot(total_time_steps-lower_value:total_time_steps+upper_value,y,'LineStyle','none','Marker','o',Color=color_order(mean_env_pair_indx));
            hold on
            plot(total_time_steps*ones(21,1),-10:10,'--k','LineWidth',2)
        end

        ax = gca;
        ax.FontSize = 16;
        xlim([analysis_time(1)-lower_value analysis_time(1)+upper_value])
        ax.XTick = analysis_time(1)-lower_value:24:analysis_time(1)+upper_value;
        ax.XTickLabel = (-lower_value:24:upper_value)/24;
        if(out_type == 1)
            ax.YTick = 0:0.2:1;
            ylim([0 1])
            ylabel('P(S_n = S^h)')
        else
            ylabel('Growth rate (h^{-1})')
            ylim([ min(cost) max(reward)])
            ax.YTick = min([reward cost]):(max([reward cost])-min([reward cost]))/4:max([reward cost]);

        end
        xlabel('time (days)');
        grid on
    end
end

%% dynamics of cell counts one-step shift in the environment and single clonal population 
%%% Numerical %%%

% num_rep= 1000;
% c0_indx = 1;
% x_num = zeros(sum(analysis_time),size(mean_env_pair,1),num_rep);
% x_num(1,:,:) = 10^4;
% 
% for mean_env_pair_indx = 1:size(mean_env_pair,1)
%     env_mat = reshape(env_time_series(:,c0_indx,mean_env_pair_indx,:),total_time_steps,num_rep);
%     state_mat = reshape(num_cell_state_time_series(:,c0_indx,mean_env_pair_indx,:),total_time_steps,num_rep);
% 
%     for rep_indx = 1:num_rep
%         for n = 1:sum(analysis_time)-1
%             % determining fitness of the clones for the current time step
%             switch env_mat(n,rep_indx)
%                 case 0
%                     if(state_mat(n,rep_indx) == 0)
%                         f = reward(2);
%                     else
%                         f = cost(1);
%                     end
%                 case 1
% 
%                     if(state_mat(n,rep_indx) == 0)
%                         f = cost(2);
%                     else
%                         f = reward(1);
%                     end
%             end
%             if(f >= 0)
%                 num_div = poissrnd(f*x_num(n,mean_env_pair_indx,rep_indx));
%                 num_death = 0;
%             else
%                 num_div = 0;
%                 num_death = poissrnd(-f*x_num(n,mean_env_pair_indx,rep_indx));
%             end
% 
%             x_num(n+1,mean_env_pair_indx,rep_indx) = x_num(n,mean_env_pair_indx,rep_indx) + num_div - num_death;
%             if(x_num(n+1,mean_env_pair_indx,rep_indx) < 0)
%                 x_num(n+1,mean_env_pair_indx,rep_indx) = 0;
%             end
%         end
%     end
% end


%%%% Analytical %%%%

x = zeros(total_time_steps,size(mean_env_pair,1));
x(1,:) = 10^4;
for mean_env_pair_indx = 1:size(mean_env_pair,1)
    for n = 1:total_time_steps-1
        % determining fitness of the clones for the current time step
        if(n<=24)
            f = reward(2);
        else
            f = benefit_time_series(n-analysis_time(1)+1,c0_indx,mean_env_pair_indx) ;
        end
        if(f >= 0)
            num_div = f*x(n,mean_env_pair_indx);
            num_death = 0;
        else
            num_div = 0;
            num_death = -f*x(n,mean_env_pair_indx);
        end

        x(n+1,mean_env_pair_indx) = x(n,mean_env_pair_indx) + num_div - num_death;
        if(x(n+1,mean_env_pair_indx) < 0)
            x(n+1,mean_env_pair_indx) = 0;
        end
    end
end

mean_x_count = x;

%%% Plotting data
color_order = ["#7E2F8E","#D95319","#EDB120","#77AC30","#0072BD"];
upper_value =  24*14;
lower_value = 24;
figure('Position',[597 568 430 236])
for mean_env_pair_indx = [size(mean_env_pair,1):-1:1]
    hold on
    plot(analysis_time(1)-lower_value:analysis_time(1)+upper_value,mean_x_count(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1,mean_env_pair_indx),'LineWidth',2,'Color',color_order(mean_env_pair_indx))
end
plot(ones(1,10)*analysis_time(1)+1,linspace(min(mean_x_count(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1,:),[],"all"),max(mean_x_count(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1,:),[],"all"),10),'--k','LineWidth',2);

xlabel('time')
ylabel('Mean cell counts')
ax = gca;
ax.FontSize = 16;
grid on
yscale log
ax.YTick = 10.^(0:11); %[10^2 10^3 10^4 10^5 10^6 10^7 10^8];
ax = gca;
xlim([analysis_time(1)-lower_value analysis_time(1)+upper_value])
ax.XTick = analysis_time(1)-lower_value:24:analysis_time(1)+upper_value;
ax.XTickLabel = (-lower_value:24:upper_value)/24;
ylim([min(mean_x_count(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1,:),[],"all") max(mean_x_count(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1,:),[],"all")])
xlabel('time (days)');
grid on
%% scatter plot of benefit in Sl state and mean time of state switching %%

sub_opt_benefit = zeros(1,size(mean_env_pair,1));
mean_adaptation_time = zeros(1,size(mean_env_pair,1));
mean_adaptation_time_num = zeros(1,size(mean_env_pair,1));

c0_indx = 1;
figure('Position',[680 597 327 281])
for mean_env_pair_indx = 1:size(mean_env_pair,1)
    temp = adaptation_time(c0_indx,mean_env_pair_indx);
    temp = temp(~isnan(temp));
    mean_adaptation_time_num(mean_env_pair_indx) = mean(temp);
    mean_adaptation_time(mean_env_pair_indx) = adaptation_time(c0_indx,mean_env_pair_indx);
    sub_opt_benefit(mean_env_pair_indx) = reward(2)*(1-mean_env_pair(mean_env_pair_indx,2))+cost(2)*(mean_env_pair(mean_env_pair_indx,2));
end

scatter(mean_adaptation_time,sub_opt_benefit,25,mean_env_pair(:,2),"filled");
hold on
color_order = ["#7E2F8E", "#D95319","#EDB120","#77AC30","#0072BD"];
mean_env_pair_plot = [0.2 0.4 0.6 0.8];
for itr = 1:length(mean_env_pair_plot)
    indx = find(round(mean_env_pair(:,2),2) == round(mean_env_pair_plot(itr),2),1);
    plot(mean_adaptation_time(indx),sub_opt_benefit(indx),'Marker','o','MarkerSize',10,'MarkerEdgeColor',color_order(itr),'LineWidth',2)
end

c = colorbar;
c.Limits = [min(mean_env_pair(:,2)) 1];
c.Ticks = min(mean_env_pair(:,2)):0.1:1;
ax = gca;
xlabel('Mean time to S^l \rightarrow S^h transition (days)')
ylabel('Benefit from S^l state')
% ax.XTick = 0:20:100;
ax.FontSize = 16;
grid on
% xscale log
ax.XTick = [0:24:max(mean_adaptation_time)+24];
ax.XTickLabel = [0:24:max(mean_adaptation_time)+24]/24;
%% Two step change in the environment and single clonal population 
num_rep = 1000;
intermediate_env = [0 0.4 0.6 0.8];

% Analytical and Numerical Calculations
c0_set = 0;%linspace(0,min(alpha,beta)-min(alpha,beta)/2,3);
num_benefit_centered_moments = zeros(4,length(c0_set),length(c0_set),length(intermediate_env),num_rep); % first two indices corresponds to the number of p and model values
num_benefit_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),length(intermediate_env),num_rep); % first two indices corresponds to the number of p and model values
num_cell_state_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),length(intermediate_env),num_rep); % first two indices corresponds to the number of p and model values
env_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),length(intermediate_env),num_rep); % first two indices corresponds to the number of p and model values

intermediate_env_length = length(intermediate_env);
end_env_pair = [0 1];

% Analytical calculation of priming duration
adaptation_time = zeros(length(c0_set),intermediate_env_length); % first two indices corresponds to the number of p and model values
c0_length = length(c0_set);
for c0_indx = 1:c0_length
    q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
    q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
    for intermediate_env_indx = 2:intermediate_env_length % starting from index 2 to avoid calulcation for p = 0
        [centered_moments,~,prob_state_1_start_0,prob_Sn_Sh_Sn_1_Sh_start_0,prob_Sn_Sl_Sn_1_Sl_start_0] = one_step_env_switch_mean_adaptation_time(N,q01,q10,end_env_pair(1),intermediate_env(intermediate_env_indx),analysis_time(2));
        adaptation_time(c0_indx,intermediate_env_indx) = centered_moments(1,1);
    end
end

%%%% dynamics of per-capita growth and cell fraction in the population 
benefit_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),intermediate_env_length); % first two indices corresponds to the number of p and model values
cell_state_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),intermediate_env_length); % first two indices corresponds to the number of p and model values

for c0_indx = 1:c0_length
    for intermediate_env_indx = 1:intermediate_env_length
        for prime_indx = 1:c0_length
            q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
            q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

            primining_duration = ceil((adaptation_time(prime_indx,intermediate_env_indx)));
            % analytical calculation of change in cell state and benefit for two step
            % change in the environment
            [prob_state_1_start_0,prob_Sn_Sh_Sn_1_Sh_start_0,prob_Sn_Sl_Sn_1_Sl_start_0] = two_step_env_switch_cell_state_dyn(N,q01,q10,end_env_pair,intermediate_env(intermediate_env_indx),analysis_time(2),primining_duration);
            % benefit_time_series(:,c0_indx,prime_indx,intermediate_env_indx) = benefit_centered_moments_cal_two_step_env_change(end_env_pair,intermediate_env(intermediate_env_indx),reward, cost, prob_Sn_Sh_Sn_1_Sh_start_0, prob_state_1_start_0,c0_set(c0_indx),primining_duration);
            % cell_state_time_series(:,c0_indx,prime_indx,intermediate_env_indx) = prob_state_1_start_0;
            
            P_s = zeros(2,2,length(prob_Sn_Sl_Sn_1_Sl_start_0));
            P_s(1,1,:) = prob_Sn_Sl_Sn_1_Sl_start_0;
            P_s(1,2,:) = 1-P_s(1,1,:);
            P_s(2,2,:) = prob_Sn_Sh_Sn_1_Sh_start_0;
            P_s(2,1,:) = 1-P_s(2,2,:);

            % Dynamics of population cell fraction and cellular benefit
            x = zeros(total_time_steps,3);
            fd = zeros(total_time_steps,1);
            x([1 2],:) = repmat([1 0 0],2,1);
            benefit_time_series(1,c0_indx,prime_indx,intermediate_env_indx) = reward(2);

            for n = 2:analysis_time(2)-1
                if(n<=primining_duration+1 && intermediate_env(intermediate_env_indx)~=0)
                    p = intermediate_env(intermediate_env_indx);
                else
                    p = end_env_pair(2);
                end

                if((1-prob_state_1_start_0(n))~=0)
                    S_l_fit = reward(2)*(1-p) +cost(2)*(p)-prob_state_1_start_0(n-1)*P_s(2,1,n)*c0_set(c0_indx)/(1-prob_state_1_start_0(n));
                else
                    S_l_fit = reward(2)*(1-p) +cost(2)*(p);
                end

                if((prob_state_1_start_0(n))~=0)
                    S_h_fit = reward(1)*(p)+cost(1)*(1-p)-(1-prob_state_1_start_0(n-1))*P_s(1,2,n)*c0_set(c0_indx)/prob_state_1_start_0(n);
                else
                    S_h_fit = reward(1)*(p)+cost(1)*(1-p);
                end

                transiton_matrix = zeros(2+1);
                growth_matrix = zeros(2+1);

                transiton_matrix(1:2, 3) = [-S_l_fit*(1-double(S_l_fit>0)); -S_h_fit*(1-double(S_h_fit>0))];
                transiton_matrix(1:2,1:2) = (1-transiton_matrix(1:2,3)).*[P_s(1,1,n) P_s(1,2,n); P_s(2,1,n) P_s(2,2,n)];
                growth_matrix(1:2,1:2) = [S_l_fit*double(S_l_fit>0) 0; 0 S_h_fit*double(S_h_fit>0)];

                transiton_matrix(3,:) = [0 0 1];

                x(n+1,:) = x(n,:)*(transiton_matrix+growth_matrix);
                fd(n) = x(n,:)*sum(transiton_matrix+growth_matrix,2);
                x(n+1,:) = x(n+1,:)/fd(n);

                benefit_time_series(n,c0_indx,prime_indx,intermediate_env_indx) = S_l_fit*(x(n,1)./sum(x(n,[1 2]),2)) + S_h_fit*(x(n,2)./sum(x(n,[1 2]),2));
            end
            cell_state_time_series(:,c0_indx,prime_indx,intermediate_env_indx) = x(:,2)./sum(x(:,[1 2]),2);

            % % Numerical Calculations
            % for rep_indx = 1:num_rep
            % rng(rep_indx)
            %     disp(['running analysis for rep num ' num2str(rep_indx) ' c0 ' num2str(c0_indx) ' intermediate_env ' num2str(intermediate_env_indx) ' prime_dur ' num2str(prime_indx) ])
            % 
            %     p = end_env_pair(1);
            %     p2 = min(1,(1-one_step_auto_corr)*(1-p));
            %     p1 = min(1,(1-one_step_auto_corr)*(p));
            % 
            %     env_1 = correlated_env(p1,p2,total_time_steps+N); % the state of the environment along each axis and timesteps
            % 
            %     p = intermediate_env(intermediate_env_indx);
            %     p2 = min(1,(1-one_step_auto_corr)*(1-p));
            %     p1 = min(1,(1-one_step_auto_corr)*(p));
            %     if(primining_duration ~= 0 )
            %         env_2 = correlated_env(p1,p2,primining_duration); % the state of the environment along each axis and timesteps
            %     else
            %         env_2 = [];
            %     end
            %     p = end_env_pair(2);
            %     p2 = min(1,(1-one_step_auto_corr)*(1-p));
            %     p1 = min(1,(1-one_step_auto_corr)*(p));
            % 
            %     env_3 = correlated_env(p1,p2,total_time_steps-primining_duration); % the state of the environment along each axis and timesteps
            % 
            %     env = [env_1 env_2 env_3]';
            % 
            %     temp_cell_state = zeros(total_time_steps,num_rep);
            %     temp_pi_est = zeros(total_time_steps,num_rep);
            %     % disp(['running numerical analyses for model ' num2str(model_indx) ' N ' num2str(N) ' c0 ' num2str(c0_set(c0_indx)) ' one step autocorr ' num2str(one_step_auto_corr(corr_env_indx)) ' rep ' num2str(rep_indx)])
            % 
            %     [cell_state_with_time,benefit,avg_benefit,avg_benefit_with_time,~] = Phenotypic_adaptation_models(reward,cost,2*total_time_steps, env, N, p, c0_set(c0_indx));
            %     temp_auto_corr = autocorr(benefit);
            % 
            %     num_benefit_centered_moments(:,c0_indx,prime_indx,intermediate_env_indx,rep_indx) = [avg_benefit; var(benefit); skewness(benefit); temp_auto_corr(2)]; % saving lag 1 autocorrelation
            %     num_benefit_time_series(:,c0_indx,prime_indx,intermediate_env_indx,rep_indx) = benefit;
            %     num_cell_state_time_series(:,c0_indx,prime_indx,intermediate_env_indx,rep_indx) = cell_state_with_time;
            %     env_time_series(:,c0_indx,prime_indx,intermediate_env_indx,rep_indx) = env(N+1:end);
            % 
            % end
        end
    end
end

%% Dynamics of per-capita growth rate and changes in the phenotypic structure of clonal population
lower_value = 24;
upper_value = 10*24;

% color_order = ["#0072BD","#7E2F8E", "#D95319","#EDB120","#77AC30"];
color_order = ["#0072BD", "#D95319","#EDB120","#77AC30","#7E2F8E"];
line_style = {'-','--'};

for out_type = 1:2 % 1:cell state, 2:avg benefit
    figure('Position',[597 568 430 236])
    for intermediate_env_indx = [1 length(intermediate_env):-1:2]
        for c0_indx = 1
            for prime_indx = 1
                primining_duration = ceil(adaptation_time(prime_indx,intermediate_env_indx));

                if(out_type == 1)
                    % data = reshape(num_cell_state_time_series(:,c0_indx,prime_indx,intermediate_env_indx,:),2*total_time_steps,num_rep);
                    y_analy = [zeros(analysis_time(1),1); cell_state_time_series(2:end,c0_indx,prime_indx,intermediate_env_indx)];
                else
                    % data = reshape(num_benefit_time_series(:,c0_indx,prime_indx,intermediate_env_indx,:),2*total_time_steps,num_rep);
                    y_analy = [ones(analysis_time(1),1)*reward(2); benefit_time_series(2:end,c0_indx,prime_indx,intermediate_env_indx)];
                end
                plot(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value,y_analy(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value),'LineWidth',2,'LineStyle','-','Color',color_order(intermediate_env_indx));
                hold on
                plot(ones(1,10)*analysis_time(1)+primining_duration+1,linspace(-1,1,10),'--','Color',color_order(intermediate_env_indx),'LineWidth',2);
            end
        end
        xlim([analysis_time(1)-lower_value analysis_time(1)+upper_value])
        ax = gca;
        ax.FontSize = 16;
        ax.XTick = analysis_time(1)-lower_value+1:24:analysis_time(1)+upper_value;
        ax.XTickLabel = (-lower_value:24:upper_value)/24;

        if(out_type == 1)
            ax.YTick = 0:0.2:1;
            ylim([0 1])
            ylabel('P(S_n = S^h)')
        else
            ax.YTick = min([reward cost]):(max([reward cost])-min([reward cost]))/4:max([reward cost]);
            ylim([min([reward cost]) max([reward cost])])
            ylabel('Growth rate (h^{-1})')
        end
        xlabel('time (days)');
        grid on
    end
end

%% dynamics of cell counts

% % Numerical
% c0_indx = 1;
% prime_indx = 1;
% 
% analysis_time = 2*total_time_steps;
% x_num = zeros(sum(analysis_time),size(mean_env_pair,1),num_rep);
% x_num(1,:,:) = 10^4;
% 
% for intermediate_env_indx = 1:length(intermediate_env)
%     env_mat = reshape(env_time_series(:,c0_indx,prime_indx,intermediate_env_indx,:),2*total_time_steps,num_rep);
%     state_mat = reshape(num_cell_state_time_series(:,c0_indx,prime_indx,intermediate_env_indx,:),2*total_time_steps,num_rep);
% 
%     for rep_indx = 1:num_rep
%         for n = 1:sum(analysis_time)-1
%             % determining fitness of the clones for the current time step
%             switch env_mat(n,rep_indx)
%                 case 0
%                     if(state_mat(n,rep_indx) == 0)
%                         f = reward(2);
%                     else
%                         f = cost(1);
%                     end
%                 case 1
% 
%                     if(state_mat(n,rep_indx) == 0)
%                         f = cost(2);
%                     else
%                         f = reward(1);
%                     end
%             end
%             if(f >= 0)
%                 num_div = poissrnd(f*x_num(n,intermediate_env_indx,rep_indx));
%                 num_death = 0;
%             else
%                 num_div = 0;
%                 num_death = poissrnd(-f*x_num(n,intermediate_env_indx,rep_indx));
%             end
% 
%             x_num(n+1,intermediate_env_indx,rep_indx) = x_num(n,intermediate_env_indx,rep_indx) + num_div - num_death;
%             if(x_num(n+1,intermediate_env_indx,rep_indx) < 0)
%                 x_num(n+1,intermediate_env_indx,rep_indx) = 0;
%             end
%         end
%     end
% end

% analytical dynamics of cell counts of single clonal population

x = zeros(sum(analysis_time),size(mean_env_pair,1));
x(1,:) = 10^4;
for intermediate_env_indx = 1:length(intermediate_env)
    for n = 1:sum(analysis_time)-1
        % determining fitness of the clones for the current time step
        if(n<=analysis_time(1))
            f = reward(2);
        else
            f = benefit_time_series(n-analysis_time(1)+1,c0_indx,prime_indx,intermediate_env_indx) ;
        end
        if(f >= 0)
            num_div = f*x(n,intermediate_env_indx);
            num_death = 0;
        else
            num_div = 0;
            num_death = -f*x(n,intermediate_env_indx);
        end

        x(n+1,intermediate_env_indx) = x(n,intermediate_env_indx) + num_div - num_death;
        if(x(n+1,intermediate_env_indx) < 0)
            x(n+1,intermediate_env_indx) = 0;
        end
    end
end

% mean_x_count = mean(x_num,3);
mean_x_count = x;

%%% Plotting data
marker_type = {'o','square','d','>','<',};
color_order = ["#0072BD", "#D95319","#EDB120","#77AC30","#7E2F8E"];
lower_value = 24;
upper_value = 14*24;
figure('Position',[597 568 430 236])
for intermediate_env_indx = [1 length(intermediate_env):-1:2]
    hold on
    plot(analysis_time(1)-lower_value:analysis_time(1)+upper_value,mean_x_count(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1,intermediate_env_indx),'LineWidth',2,'Color',color_order(intermediate_env_indx))
    plot(ones(1,100)*analysis_time(1)+adaptation_time(1,intermediate_env_indx)+1,linspace(min(mean_x_count,[],"all"),max(mean_x_count,[],"all"),100),'LineWidth',2,'LineStyle','--','Color',color_order(intermediate_env_indx));
end
xlabel('time')
ylabel('Mean cell counts')
ax = gca;
ax.FontSize = 16;
grid on
yscale log
ylim([min(mean_x_count(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1,1:4),[],"all") max(mean_x_count(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1,1:4),[],"all")])
ax.YTick = 10.^(0:11); %[10^2 10^3 10^4 10^5 10^6 10^7 10^8];
ax = gca;
xlim([analysis_time(1)-lower_value analysis_time(1)+upper_value])
ax.XTick = [analysis_time(1)-lower_value:24:analysis_time(1)+upper_value];
ax.XTickLabel = [-lower_value:24:upper_value]/24;
xlabel('time (days)');
grid on

%% mannual variation of priming duration
analysis_time = [24 24*28];
total_time_steps = sum(analysis_time);

intermediate_env = [0 0.4 0.6 0.8];
c0_set = linspace(0,min(alpha,beta)/2,3);
N = 20:20:60;
intermediate_env_length = length(intermediate_env);
end_env_pair = [0 1];

priming_duration_set = 24:24:24*14;

benefit_time_series = zeros(analysis_time(2),length(c0_set),length(priming_duration_set),intermediate_env_length); % first two indices corresponds to the number of p and model values
cell_state_time_series = zeros(analysis_time(2),length(c0_set),length(priming_duration_set),intermediate_env_length); % first two indices corresponds to the number of p and model values
figure('Units','normalized','Position',[0 0 1 1])

for N_indx = 1:length(N)
    for c0_indx = 1:length(c0_set)
        for intermediate_env_indx = 1:intermediate_env_length
            for prime_indx = 1:length(priming_duration_set)
                q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
                q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

                if(intermediate_env_indx == 1)
                    primining_duration = 0;
                else
                    primining_duration = priming_duration_set(prime_indx);
                end
                % analytical calculation of change in cell state fraction and per-capita growth for two step
                % change in the environment
                [prob_state_1_start_0,prob_Sn_Sh_Sn_1_Sh_start_0,prob_Sn_Sl_Sn_1_Sl_start_0] = two_step_env_switch_cell_state_dyn(N(N_indx),q01,q10,end_env_pair,intermediate_env(intermediate_env_indx),analysis_time(2),primining_duration);
                P_s = zeros(2,2,length(prob_Sn_Sl_Sn_1_Sl_start_0));
                P_s(1,1,:) = prob_Sn_Sl_Sn_1_Sl_start_0;
                P_s(1,2,:) = 1-P_s(1,1,:);
                P_s(2,2,:) = prob_Sn_Sh_Sn_1_Sh_start_0;
                P_s(2,1,:) = 1-P_s(2,2,:);

                % Dynamics of population cell fraction and cellular benefit
                x = zeros(analysis_time(2),3);
                fd = zeros(analysis_time(2),1);
                x([1 2],:) = repmat([1 0 0],2,1);
                benefit_time_series(1,c0_indx,prime_indx,intermediate_env_indx) = reward(2);

                for n = 2:analysis_time(2)-1
                    if(n<=primining_duration+1 && intermediate_env(intermediate_env_indx)~=0)
                        p = intermediate_env(intermediate_env_indx);
                    else
                        p = end_env_pair(2);
                    end

                    if((1-prob_state_1_start_0(n))~=0)
                        S_l_fit = reward(2)*(1-p) +cost(2)*(p)-prob_state_1_start_0(n-1)*P_s(2,1,n)*c0_set(c0_indx)/(1-prob_state_1_start_0(n));
                    else
                        S_l_fit = reward(2)*(1-p) +cost(2)*(p);
                    end

                    if((prob_state_1_start_0(n))~=0)
                        S_h_fit = reward(1)*(p)+cost(1)*(1-p)-(1-prob_state_1_start_0(n-1))*P_s(1,2,n)*c0_set(c0_indx)/prob_state_1_start_0(n);
                    else
                        S_h_fit = reward(1)*(p)+cost(1)*(1-p);
                    end

                    transiton_matrix = zeros(2+1);
                    growth_matrix = zeros(2+1);

                    transiton_matrix(1:2, 3) = [-S_l_fit*(1-double(S_l_fit>0)); -S_h_fit*(1-double(S_h_fit>0))];
                    transiton_matrix(1:2,1:2) = (1-transiton_matrix(1:2,3)).*[P_s(1,1,n) P_s(1,2,n); P_s(2,1,n) P_s(2,2,n)];
                    growth_matrix(1:2,1:2) = [S_l_fit*double(S_l_fit>0) 0; 0 S_h_fit*double(S_h_fit>0)];

                    transiton_matrix(3,:) = [0 0 1];

                    x(n+1,:) = x(n,:)*(transiton_matrix+growth_matrix);
                    fd(n) = x(n,:)*sum(transiton_matrix+growth_matrix,2);
                    x(n+1,:) = x(n+1,:)/fd(n);

                    benefit_time_series(n,c0_indx,prime_indx,intermediate_env_indx) = S_l_fit*(x(n,1)./sum(x(n,[1 2]),2)) + S_h_fit*(x(n,2)./sum(x(n,[1 2]),2));
                end
                cell_state_time_series(:,c0_indx,prime_indx,intermediate_env_indx) = x(:,2)./sum(x(:,[1 2]),2);
            end
        end

        mean_x_count_end = zeros(length(priming_duration_set),intermediate_env_length);

        % dynamics of total cell counts
        for intermediate_env_indx = 1:length(intermediate_env)
            for prime_indx = 1:length(priming_duration_set)
                x = zeros(sum(analysis_time),1);
                x(1,:) = 10^4;
                for n = 1:sum(analysis_time)-1
                    % determining fitness of the clones for the current time step
                    if(n<=analysis_time(1))
                        f = reward(2);
                    else
                        f = benefit_time_series(n-analysis_time(1)+1,c0_indx,prime_indx,intermediate_env_indx) ;
                    end
                    if(f >= 0)
                        num_div = f*x(n);
                        num_death = 0;
                    else
                        num_div = 0;
                        num_death = -f*x(n);
                    end

                    x(n+1) = x(n) + num_div - num_death;
                    if(x(n+1) < 0)
                        x(n+1) = 0;
                    end
                end
                mean_x_count_end(prime_indx,intermediate_env_indx) = x(n+1);
            end

        end

        % color_order = ["#0072BD", "#7E2F8E"	, "#D95319","#EDB120","#77AC30"];
        color_order = ["#0072BD", "#D95319","#EDB120","#77AC30"];

        subplot(3,3,3*(N_indx-1)+c0_indx)
        for intermediate_env_indx = 2:4
            plot(priming_duration_set,mean_x_count_end(:,intermediate_env_indx)./mean_x_count_end(:,1),'LineWidth',2,'Color', color_order(intermediate_env_indx));
            hold on
        end
        plot(priming_duration_set,ones(1,length(priming_duration_set)),'LineWidth',2,'Color',color_order(1),'LineStyle','--')
        xlabel('time')
        ylabel('Cell counts fold change')
        ax = gca;
        ax.FontSize = 16;
        grid on
        % yscale log
        ax = gca;
        xlim([priming_duration_set(1) priming_duration_set(end)])
        ax.XTick = [priming_duration_set(1):24:priming_duration_set(end)];
        ax.XTickLabel = [priming_duration_set(1):24:priming_duration_set(end)]/24;
        xlabel('priming duration (days)');
        % ax.YTick = 10.^(-1:3);
        grid on
    end
end
%% Two step change in the environment with multi clonal population
num_rep = 1000;
intermediate_env = [0 0.4 0.6 0.8];
analysis_time = [24 24*7];
total_time_steps = sum(analysis_time);
N = 40;

% Analytical and Numerical Calculations

% c0_set = [0.0248];
c0_set = [0.0248    0.0495];

num_benefit_centered_moments = zeros(4,length(c0_set),length(c0_set),length(intermediate_env),num_rep); % first two indices corresponds to the number of p and model values
num_benefit_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),length(intermediate_env),num_rep); % first two indices corresponds to the number of p and model values
num_cell_state_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),length(intermediate_env),num_rep); % first two indices corresponds to the number of p and model values
env_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),length(intermediate_env),num_rep); % first two indices corresponds to the number of p and model values

intermediate_env_length = length(intermediate_env);
end_env_pair = [0 1];

% Analytical calculation of priming duration
adaptation_time = 0;
primining_duration = analysis_time(1);
benefit_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),intermediate_env_length);
pop_benefit_time_series = zeros(total_time_steps,length(c0_set),intermediate_env_length);

benefit_time_series_individual_clones = zeros(total_time_steps,length(c0_set),length(c0_set),intermediate_env_length); % first two indices corresponds to the number of p and model values
cell_state_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),intermediate_env_length); % first two indices corresponds to the number of p and model values
clonal_fraction_time_series = zeros(total_time_steps,length(c0_set),length(c0_set),intermediate_env_length);

for intermediate_env_indx = 1%:2%intermediate_env_length
    for prime_indx = 1%:c0_length
        prob_state_1_start_0 = zeros(total_time_steps+1,length(c0_set));
        P_s = zeros(2,2,total_time_steps+2,length(c0_set));
        for c0_indx = 1:length(c0_set)
            q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
            q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
            % analytical calculation of change in cell state and benefit for two step
            % change in the environment
            [prob_state_1_start_0(:,c0_indx),prob_Sn_Sh_Sn_1_Sh_start_0,prob_Sn_Sl_Sn_1_Sl_start_0] = two_step_env_switch_cell_state_dyn(N,q01,q10,end_env_pair,intermediate_env(intermediate_env_indx),total_time_steps,primining_duration);
            % benefit_time_series(:,c0_indx,prime_indx,intermediate_env_indx) = benefit_centered_moments_cal_two_step_env_change(end_env_pair,intermediate_env(intermediate_env_indx),reward, cost, prob_Sn_Sh_Sn_1_Sh_start_0, prob_state_1_start_0,c0_set(c0_indx),primining_duration);
            % cell_state_time_series(:,c0_indx,prime_indx,intermediate_env_indx) = prob_state_1_start_0;
            P_s(1,1,:,c0_indx) = prob_Sn_Sl_Sn_1_Sl_start_0;
            P_s(1,2,:,c0_indx) = 1-P_s(1,1,:,c0_indx);
            P_s(2,2,:,c0_indx) = prob_Sn_Sh_Sn_1_Sh_start_0;
            P_s(2,1,:,c0_indx) = 1-P_s(2,2,:,c0_indx);
        end

        % Dynamics of population cell fraction and cellular benefit
        x = zeros(total_time_steps,2*length(c0_set)+1);
        fd = zeros(total_time_steps,1);
        
        x([1 2],:) = [repmat([0.004 0 0.9996 0],2,1) zeros(2,1)];
        % x([1 2],:) = [repmat([1 0],2,1) zeros(2,1)];


        clonal_fraction_time_series(1,:,prime_indx,intermediate_env_indx) = x(1,1:2:2*length(c0_set));
        benefit_time_series(1,:,prime_indx,intermediate_env_indx) = reward(2);
        benefit_time_series_individual_clones(1,:,prime_indx,intermediate_env_indx) = reward(2);

        for n = 2:total_time_steps-1
            if(n<=analysis_time(1)+1)
                p = intermediate_env(intermediate_env_indx);
            else
                p = end_env_pair(2);
            end
            transiton_matrix = zeros(2*length(c0_set)+1);
            growth_matrix = zeros(2*length(c0_set)+1);

            S_l_fit = zeros(1,length(c0_set));
            S_h_fit = zeros(1,length(c0_set));
            for c0_indx = 1:length(c0_set)

                if((1-prob_state_1_start_0(n,c0_indx))~=0)
                    S_l_fit(c0_indx) = reward(2)*(1-p) +cost(2)*(p)-prob_state_1_start_0(n-1,c0_indx)*P_s(2,1,n,c0_indx)*c0_set(c0_indx)/(1-prob_state_1_start_0(n,c0_indx));
                else
                    S_l_fit(c0_indx) = reward(2)*(1-p) +cost(2)*(p);
                end

                if((prob_state_1_start_0(n,c0_indx))~=0)
                    S_h_fit(c0_indx) = reward(1)*(p)+cost(1)*(1-p)-(1-prob_state_1_start_0(n-1,c0_indx))*P_s(1,2,n,c0_indx)*c0_set(c0_indx)/prob_state_1_start_0(n,c0_indx);
                else
                    S_h_fit(c0_indx) = reward(1)*(p)+cost(1)*(1-p);
                end

                % if(c0_indx == 2)
                %     S_h_fit(c0_indx) = S_l_fit(c0_indx); 
                % end

                clonal_fraction_time_series(n,c0_indx,prime_indx,intermediate_env_indx) =  sum(x(n,2*(c0_indx)-1:2*(c0_indx)))./sum(x(n,1:end-1),2);
                % cell_state_time_series(n,c0_indx,prime_indx,intermediate_env_indx) = (x(n,2*(c0_indx))./sum(x(n,2*(c0_indx)-1:2*(c0_indx))));

                benefit_time_series(n,c0_indx,prime_indx,intermediate_env_indx) = S_l_fit(c0_indx)*(x(n,2*(c0_indx)-1)./sum(x(n,2*(c0_indx)-1:2*(c0_indx)))) + S_h_fit(c0_indx)*(x(n,2*(c0_indx))./sum(x(n,2*(c0_indx)-1:2*(c0_indx))));
                pop_benefit_time_series(n,prime_indx,intermediate_env_indx) = pop_benefit_time_series(n,prime_indx,intermediate_env_indx) + S_l_fit(c0_indx)*(x(n,2*(c0_indx)-1)./sum(x(n,1:end-1))) + S_h_fit(c0_indx)*(x(n,2*(c0_indx))./sum(x(n,1:end-1)));

                transiton_matrix(2*(c0_indx)-1:2*(c0_indx), 2*(length(c0_set))+1) = [-S_l_fit(c0_indx)*(1-double(S_l_fit(c0_indx)>0)); -S_h_fit(c0_indx)*(1-double(S_h_fit(c0_indx)>0))];                
                transiton_matrix(2*(c0_indx)-1:2*(c0_indx),2*(c0_indx)-1:2*(c0_indx)) = (1-transiton_matrix(2*(c0_indx)-1:2*(c0_indx), 2*(length(c0_set))+1)).*[P_s(1,1,n,c0_indx) P_s(1,2,n,c0_indx); P_s(2,1,n,c0_indx) P_s(2,2,n,c0_indx)];
                growth_matrix(2*(c0_indx)-1:2*(c0_indx),2*(c0_indx)-1:2*(c0_indx)) = [S_l_fit(c0_indx)*double(S_l_fit(c0_indx)>0) 0; 0 S_h_fit(c0_indx)*double(S_h_fit(c0_indx)>0)];
            end
                transiton_matrix(end,:) = [zeros(1,2*length(c0_set)) 1];

                x(n+1,:) = x(n,:)*(transiton_matrix+growth_matrix);
                fd(n) = x(n,:)*sum(transiton_matrix+growth_matrix,2);
                x(n+1,:) = x(n+1,:)/fd(n);

        end
        cell_state_time_series(:,:,prime_indx,intermediate_env_indx) = x(:,2:2:2*length(c0_set))./sum(x(:,1:end-1),2);

        % % Numerical Calculations
        % for rep_indx = 1:num_rep
        % rng(rep_indx)
        %     disp(['running analysis for rep num ' num2str(rep_indx) ' c0 ' num2str(c0_indx) ' intermediate_env ' num2str(intermediate_env_indx) ' prime_dur ' num2str(prime_indx) ])
        %
        %     p = end_env_pair(1);
        %     p2 = min(1,(1-one_step_auto_corr)*(1-p));
        %     p1 = min(1,(1-one_step_auto_corr)*(p));
        %
        %     env_1 = correlated_env(p1,p2,total_time_steps+N); % the state of the environment along each axis and timesteps
        %
        %     p = intermediate_env(intermediate_env_indx);
        %     p2 = min(1,(1-one_step_auto_corr)*(1-p));
        %     p1 = min(1,(1-one_step_auto_corr)*(p));
        %     if(primining_duration ~= 0 )
        %         env_2 = correlated_env(p1,p2,primining_duration); % the state of the environment along each axis and timesteps
        %     else
        %         env_2 = [];
        %     end
        %     p = end_env_pair(2);
        %     p2 = min(1,(1-one_step_auto_corr)*(1-p));
        %     p1 = min(1,(1-one_step_auto_corr)*(p));
        %
        %     env_3 = correlated_env(p1,p2,total_time_steps-primining_duration); % the state of the environment along each axis and timesteps
        %
        %     env = [env_1 env_2 env_3]';
        %
        %     temp_cell_state = zeros(total_time_steps,num_rep);
        %     temp_pi_est = zeros(total_time_steps,num_rep);
        %     % disp(['running numerical analyses for model ' num2str(model_indx) ' N ' num2str(N) ' c0 ' num2str(c0_set(c0_indx)) ' one step autocorr ' num2str(one_step_auto_corr(corr_env_indx)) ' rep ' num2str(rep_indx)])
        %
        %     [cell_state_with_time,benefit,avg_benefit,avg_benefit_with_time,~] = Phenotypic_adaptation_models(reward,cost,2*total_time_steps, env, N, p, c0_set(c0_indx));
        %     temp_auto_corr = autocorr(benefit);
        %
        %     num_benefit_centered_moments(:,c0_indx,prime_indx,intermediate_env_indx,rep_indx) = [avg_benefit; var(benefit); skewness(benefit); temp_auto_corr(2)]; % saving lag 1 autocorrelation
        %     num_benefit_time_series(:,c0_indx,prime_indx,intermediate_env_indx,rep_indx) = benefit;
        %     num_cell_state_time_series(:,c0_indx,prime_indx,intermediate_env_indx,rep_indx) = cell_state_with_time;
        %     env_time_series(:,c0_indx,prime_indx,intermediate_env_indx,rep_indx) = env(N+1:end);
        %
        % end
    end
end
% figure
% subplot(4,1,1)
% temp = clonal_fraction_time_series(:,:,1,1);
% plot(1:total_time_steps,temp(1:total_time_steps,:),'LineWidth',2');
% ylim([0.2 0.45])
% intermediate_env_indx = 2;
% for prime_indx = 1:length(c0_set)
% subplot(4,1,prime_indx+1)
% temp = clonal_fraction_time_series(:,:,prime_indx,intermediate_env_indx);
% plot(1:total_time_steps,temp(1:total_time_steps,:),'LineWidth',2);
% ylim([0.2 0.45])
% end
%% dynamic of per-capita growth rates of individual clonal population

marker_type = {'square','d','>','<',};
max_primining_duration = ceil(max(adaptation_time,[],'all'));
lower_value = 24;
upper_value = 24*6.99;
intermediate_env_indx = 2;%:length(intermediate_env)

for out_type = 2%1:2%
    figure('Position',[1069 555 466 251])
    color_order = ["#0072BD", "#D95319","#EDB120","#77AC30","#7E2F8E"];
    tiledlayout(1,1,"TileSpacing","compact","Padding","loose");

    % tiledlayout(4,1,"TileSpacing","compact","Padding","loose");
    nexttile
    for c0_indx = 1:length(c0_set)
        if(out_type == 1)
            % data = reshape(num_cell_state_time_series(:,c0_indx,1,1,:),2*total_time_steps,num_rep);
            % y_analy = [zeros(total_time_steps,1); cell_state_time_series(2:end,c0_indx,1,1)];
            y_analy = [cell_state_time_series(1:end,c0_indx,1,1)];

        else
            % data = reshape(num_benefit_time_series(:,c0_indx,1,1,:),2*total_time_steps,num_rep);
            y_analy = benefit_time_series(1:end,c0_indx,1,1);
            y_analy(end) = y_analy(end-1); 
        end
        % y = mean(data(total_time_steps-lower_value+1:total_time_steps+upper_value,:),2);
        hold on
        plot(analysis_time(1)-lower_value:analysis_time(1)+upper_value,y_analy(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1),'LineWidth',2,'LineStyle','-','Color',color_order(c0_indx));
        % plot(total_time_steps-lower_value+1:2:total_time_steps+upper_value,y(1:2:end),'LineStyle','none','Marker','o','Color',color_order(c0_indx));

        plot(ones(1,10)*analysis_time(1)+1,linspace(-1,1,10),'--k','LineWidth',2);
    end

    % xlim([total_time_steps-lower_value total_time_steps+upper_value])
    ax = gca;
    ax.FontSize = 14;
    xlim([analysis_time(1)-lower_value analysis_time(1)+floor(upper_value/24)*24+24])
    ax.XTick = (analysis_time(1)-lower_value:24:analysis_time(1)+floor(upper_value/24)*24+24);
    ax.XTickLabel = (-lower_value:24:floor(upper_value/24)*24+24)/24;
    xlabel('time (days)')
    if(out_type == 1)
        ax.YTick = 0:0.2:1;
        ylim([0 1])
        ylabel('P(S_n = S^h)')
    else
        ax.YTick = min([reward cost]):(max([reward cost])-min([reward cost]))/4:max([reward cost]);
        ylim([min([reward cost])-0.005 max([reward cost])])
        ylabel('Growth rate (h^{-1})')
    end
    grid on
end

%% analytical dynamics of cell counts of multi-clonal population
intermediate_env_priming_indices_pair = [1 1;
    2 1;
    2 2;
    2 3];
%%% Plotting data
% color_order = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];
lower_value = 24;
upper_value = 24*6.99;
%%%%% Changes in the clone size
color_order = ["#0072BD", "#D95319","#EDB120","#7E2F8E"];

figure('Position',[1069 555 466 251])
% figure('Position',[597 75 466 729])

% tiledlayout(length(c0_set)+1,1,"TileSpacing","compact","Padding","loose");
tiledlayout(1,1,"TileSpacing","compact","Padding","loose");

for indices_pair = 1%:length(c0_set)+1
    prime_indx = intermediate_env_priming_indices_pair(indices_pair,2);
    intermediate_env_indx = intermediate_env_priming_indices_pair(indices_pair,1);
    %%%%% population dynamics %%%%%%%%%
    x = zeros(sum(analysis_time),length(c0_set));
    x(1,:) = 10^4*(clonal_fraction_time_series(1,:,prime_indx,intermediate_env_indx));

    for n = 1:sum(analysis_time)-1
        for c0_indx = 1:length(c0_set)
            % if(n<=total_time_steps)
            %     f = reward(2);
            % else
                f = (benefit_time_series(n+1,c0_indx,prime_indx,intermediate_env_indx));
            % end
            if(f >= 0)
                num_div = f*x(n,c0_indx);
                num_death = 0;
            else
                num_div = 0;
                num_death = -f*x(n,c0_indx);
            end

            x(n+1,c0_indx) = x(n,c0_indx) + num_div - num_death;
            if(x(n+1,c0_indx) < 0)
                x(n+1,c0_indx) = 0;
            end
        end
    end

    mean_x_count = x;
    x_frac = x./sum(x,2);
    mean_cell_state_frac = x_frac.*[cell_state_time_series(1:total_time_steps,:,prime_indx,intermediate_env_indx)];
    mean_x_frac = x_frac;

    x_tot = zeros(sum(analysis_time),1);
    x_tot(1,:) = 10^4;

    for n = 1:sum(analysis_time)-1

    f = pop_benefit_time_series(n+1,prime_indx,intermediate_env_indx);

            if(f >= 0)
                num_div = f*x_tot(n);
                num_death = 0;
            else
                num_div = 0;
                num_death = -f*x_tot(n);
            end

            x_tot(n+1) = x_tot(n) + num_div - num_death;
            if(x_tot(n+1) < 0)
                x_tot(n+1) = 0;
            end
        
    end

    mean_x_tot_count = x_tot;

    
    for c0_indx = 1:length(c0_set)
        nexttile(indices_pair)

        hold on
        plot(analysis_time(1)-lower_value:analysis_time(1)+upper_value,mean_x_count(analysis_time(1)-lower_value+1:analysis_time(1)+upper_value+1,c0_indx),'LineWidth',2,'Color',color_order(c0_indx))
        hold on
        if(c0_indx == length(c0_set))
            plot(ones(1,20)*analysis_time(1)+1,linspace(10^-2,1*10^8,20),'k','LineWidth',2,'LineStyle','--');
        end

        if(indices_pair == length(c0_set)+1)
                xlabel('time')
        end
        ylabel('mean cell counts')
        ax = gca;
        ax.FontSize = 14;
        grid on
        yscale log
        ax.YTick = 10.^(-1:6); %[10^2 10^3 10^4 10^5 10^6 10^7 10^8];
        % xlim([total_time_steps-lower_value total_time_steps+upper_value])
        ax = gca;
        ax.FontSize = 14;
        xlim([analysis_time(1)-lower_value analysis_time(1)+floor(upper_value/24)*24+24])
        ax.XTick = (analysis_time(1)-lower_value:24:analysis_time(1)+floor(upper_value/24)*24+24);
        ax.XTickLabel = (-lower_value:24:floor(upper_value/24)*24+24)/24;
        ylim([10^-1 10^5])
        xlabel('time (days)')

        %%%%% Changes in the clone fraction
        % nexttile((indices_pair-1)*3+1)
        % nexttile(indices_pair)
        % hold on
        % plot(total_time_steps-lower_value:total_time_steps+upper_value,mean_x_frac(total_time_steps-lower_value+1:total_time_steps+upper_value+1,c0_indx),'LineWidth',2,'Color',color_order(c0_indx))
        % plot(ones(1,10)*total_time_steps+1,linspace(-1,1,10),'--k','LineWidth',2);
        % plot(ones(1,10)*total_time_steps+adaptation_time(prime_indx,intermediate_env_indx)+1,linspace(-1,1,10),'--k','LineWidth',2);
        % 
        % if(indices_pair == length(c0_set)+1)
        %         xlabel('time')
        % end
        % ylabel('Subpop fraction')
        % ax = gca;
        % ax.FontSize = 14;
        % grid on
        % % ax.YTick = 0.25:0.05:0.4;
        % % ylim([min(mean_x_frac_analy,[],"all") max(mean_x_frac_analy,[],"all")])
        % xlim([total_time_steps-lower_value total_time_steps+upper_value])
        % ax = gca;
        % ax.FontSize = 14;
        % ax.XTick = total_time_steps-lower_value+1:20:total_time_steps+upper_value;
        % ax.XTickLabel = -lower_value:20:upper_value;
        % ylim([0 1])

        % nexttile((indices_pair-1)*3+3)
        % hold on
        % plot(total_time_steps-lower_value:total_time_steps+upper_value,mean_cell_state_frac(total_time_steps-lower_value+1:total_time_steps+upper_value+1,c0_indx),'LineWidth',2)
        % if(c0_indx == length(c0_set))
        %     plot(total_time_steps-lower_value:total_time_steps+upper_value,sum(mean_cell_state_frac(total_time_steps-lower_value+1:total_time_steps+upper_value+1,:),2),'--','LineWidth',2)
        % end
        % xlabel('time')
        % ylabel('Cell state fraction')
        % ax = gca;
        % ax.FontSize = 14;
        % grid on
        % xlim([total_time_steps-lower_value total_time_steps+upper_value])
        % ax = gca;
        % ax.FontSize = 14;
        % ax.XTick = total_time_steps-lower_value+1:10:total_time_steps+upper_value;
        % ax.XTickLabel = -lower_value:10:upper_value;
        % ax.YTick = 0:0.2:1;
        % % ylim([min(mean_x_frac_analy,[],"all") max(mean_x_frac_analy,[],"all")])
        % % plt_indx = plt_indx+1;
    end
end

%% Variation of adaptation time with increasing cost and memory size
total_time_steps = 150;

cost = [0 -1].*(1/20);
reward = [0.2 1].*(1/20);

q = (reward(2)-cost(1))/(sum(reward) - sum(cost));
alpha = reward(1)-cost(2);
beta = alpha*q./(1-q); % indifference environment (p_I)

mean_env_pair = [0 1];
N = 20:60;
c0_set = linspace(0,min(alpha,beta)-min(alpha,beta)/100,3);

adaptation_time = zeros(length(c0_set),length(N)); % first two indices corresponds to the number of p and model values

c0_length = length(c0_set);
N_length = length(N);

for c0_indx = 1:c0_length
    q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
    q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
    for N_indx = 1:N_length
        disp(['running analysis for c0 ' num2str(c0_indx) ' N ' num2str(N_indx) ])
        p = mean_env_pair(2);
        [centered_moments,~,prob_state_1_start_0,prob_Sn_Sh_Sn_1_Sh_start_0,prob_Sn_Sl_Sn_1_Sl_start_0] = one_step_env_switch_mean_adaptation_time(N(N_indx),q01,q10,mean_env_pair(1),mean_env_pair(2),total_time_steps);
        adaptation_time(c0_indx,N_indx) = centered_moments(1,1);
    end
end

figure('Position',[1069 555 466 251])
for c0_indx = 2:3
    hold on
    plot(N,adaptation_time(c0_indx,:),'LineWidth',2);
end
% yscale log
ax = gca;
ax.FontSize = 14;
grid on
xlabel('memory size (m)');
ylabel([{'Mean time to S^l \rightarrow S^h'}, {'transition (days)'}])
xlim([min(N) max(N)])
% for c0_indx = 1
ax.YTick = [0:24:max(adaptation_time,[],'all')+24];
ax.YTickLabel = [0:24:max(adaptation_time,[],'all')+24]/24;

%% Effect of memory size on advantage of quick adaptation with smaller adaptation cost
intermediate_env = 0;
% analysis_time = [24 24*14];
analysis_time = [24 24*7];
N = 20:60;
total_time_steps = sum(analysis_time);

end_env_pair = [0 1];
benefit_time_series = zeros(analysis_time(2),length(N),2);
mean_cell_count = zeros(total_time_steps,length(N),2);

for sce_indx = 1:2
    for N_indx = 1:length(N)
        switch sce_indx 
            case 1
                % c0_set = [0 0.0248    0.0495];
                c0_set = [0.0248    0.0495];
                x = zeros(analysis_time(2),2*length(c0_set)+1);
                x([1 2],:) = [repmat([0.004 0 0.9996 0],2,1) zeros(2,1)];

            case 2
                % c0_set = [0 0.0248];
                c0_set = [0.0248];

                x = zeros(analysis_time(2),2*length(c0_set)+1);
                % x([1 2],:) = [repmat([0.2 0 0.8 0],2,1) zeros(2,1)];
                x([1 2],:) = [repmat([1 0],2,1) zeros(2,1)];

        end

        primining_duration = 0;
        prob_state_1_start_0 = zeros(analysis_time(2)+1,length(c0_set));
        P_s = zeros(2,2,analysis_time(2)+2,length(c0_set));
        for c0_indx = 1:length(c0_set)
            q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
            q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
            % analytical calculation of change in cell state and benefit for two step
            % change in the environment
            [prob_state_1_start_0(:,c0_indx),prob_Sn_Sh_Sn_1_Sh_start_0,prob_Sn_Sl_Sn_1_Sl_start_0] = two_step_env_switch_cell_state_dyn(N(N_indx),q01,q10,end_env_pair,intermediate_env,analysis_time(2),primining_duration);
            P_s(1,1,:,c0_indx) = prob_Sn_Sl_Sn_1_Sl_start_0;
            P_s(1,2,:,c0_indx) = 1-P_s(1,1,:,c0_indx);
            P_s(2,2,:,c0_indx) = prob_Sn_Sh_Sn_1_Sh_start_0;
            P_s(2,1,:,c0_indx) = 1-P_s(2,2,:,c0_indx);
        end

        % Dynamics of population cell fraction and cellular benefit
        fd = zeros(analysis_time(2),1);

        benefit_time_series(1,N_indx,sce_indx) = reward(2);
        for n = 2:analysis_time(2)-1
            p = end_env_pair(2);
            transiton_matrix = zeros(2*length(c0_set)+1);
            growth_matrix = zeros(2*length(c0_set)+1);

            S_l_fit = zeros(1,length(c0_set));
            S_h_fit = zeros(1,length(c0_set));
            temp_benefit = 0;
            for c0_indx = 1:length(c0_set)

                if((1-prob_state_1_start_0(n,c0_indx))~=0)
                    S_l_fit(c0_indx) = reward(2)*(1-p) +cost(2)*(p)-prob_state_1_start_0(n-1,c0_indx)*P_s(2,1,n,c0_indx)*c0_set(c0_indx)/(1-prob_state_1_start_0(n,c0_indx));
                else
                    S_l_fit(c0_indx) = reward(2)*(1-p) +cost(2)*(p);
                end

                if((prob_state_1_start_0(n,c0_indx))~=0)
                    S_h_fit(c0_indx) = reward(1)*(p)+cost(1)*(1-p)-(1-prob_state_1_start_0(n-1,c0_indx))*P_s(1,2,n,c0_indx)*c0_set(c0_indx)/prob_state_1_start_0(n,c0_indx);
                else
                    S_h_fit(c0_indx) = reward(1)*(p)+cost(1)*(1-p);
                end

                temp_benefit = temp_benefit + (S_l_fit(c0_indx)*(x(n,2*(c0_indx)-1)/sum(x(n,1:end-1))) + S_h_fit(c0_indx)*(x(n,2*(c0_indx))/sum(x(n,1:end-1))));

                transiton_matrix(2*(c0_indx)-1:2*(c0_indx), 2*(length(c0_set))+1) = [-S_l_fit(c0_indx)*(1-double(S_l_fit(c0_indx)>0)); -S_h_fit(c0_indx)*(1-double(S_h_fit(c0_indx)>0))];
                transiton_matrix(2*(c0_indx)-1:2*(c0_indx),2*(c0_indx)-1:2*(c0_indx)) = (1-transiton_matrix(2*(c0_indx)-1:2*(c0_indx), 2*(length(c0_set))+1)).*[P_s(1,1,n,c0_indx) P_s(1,2,n,c0_indx); P_s(2,1,n,c0_indx) P_s(2,2,n,c0_indx)];
                growth_matrix(2*(c0_indx)-1:2*(c0_indx),2*(c0_indx)-1:2*(c0_indx)) = [S_l_fit(c0_indx)*double(S_l_fit(c0_indx)>0) 0; 0 S_h_fit(c0_indx)*double(S_h_fit(c0_indx)>0)];
            end
            benefit_time_series(n,N_indx,sce_indx) = temp_benefit;
            transiton_matrix(end,:) = [zeros(1,2*length(c0_set)) 1];

            x(n+1,:) = x(n,:)*(transiton_matrix+growth_matrix);
            fd(n) = x(n,:)*sum(transiton_matrix+growth_matrix,2);
            x(n+1,:) = x(n+1,:)/fd(n);
        end

        %%%%% population dynamics %%%%%%%%%
        x = zeros(sum(analysis_time),1);
        x(1) = 10^4;

        for n = 1:sum(analysis_time)-1
            if(n<=analysis_time(1))
                f = reward(2);
            else
                f = (benefit_time_series(n-analysis_time(1)+1,N_indx,sce_indx));
            end
            if(f >= 0)
                num_div = f*x(n);
                num_death = 0;
            else
                num_div = 0;
                num_death = -f*x(n);
            end

            x(n+1) = x(n) + num_div - num_death;
            if(x(n+1) < 0)
                x(n+1) = 0;
            end 
        end
    
        mean_cell_count(:,N_indx,sce_indx) = x;
    end
end

% figure('Position',[680 535 390 343])
% for N_indx = 1%:length(N)
%     hold on
%     % plot(N,mean_cell_count(end,:,2)./mean_cell_count(end,:,1),'LineWidth',2);
%     plot(1:total_time_steps,mean_cell_count(:,N_indx,2),'LineWidth',2);
%     % hold on
%     plot(1:total_time_steps,mean_cell_count(:,N_indx,1),'LineWidth',2);
% end
% yscale log
% xlim([130 300])
% figure('Position',[680 535 390 343])
% plot(N,mean_cell_count(end,:,2)./mean_cell_count(end,:,1),'LineWidth',2);
% ax = gca;
% ax.FontSize = 14;
% grid on
% xlabel('memory size (m)');
% ylabel([{'Ratio total cell counts '}, {'- enriched to mixed pop'}]);
% xlim([min(N) max(N)])

figure('Position',[1069 555 466 251])
yyaxis left
plot(N,mean_cell_count(end,:,1),'-','LineWidth',2);
hold on
plot(N,mean_cell_count(end,:,2),'-.','LineWidth',2);
yscale log
ylabel([{'End time point'},{'total cell counts'}]);
ax = gca;
ax.YTick = 10.^(2:8);
yyaxis right
plot(N,mean_cell_count(end,:,2)./mean_cell_count(end,:,1),'LineWidth',2);
ylabel(['Fold change']);
ax = gca;
ax.FontSize = 14;
grid on
xlabel('memory size (m)');
xlim([min(N) max(N)])

    % xlabel('time')
    % ylabel('Cell counts')
    % ax = gca;
    % ax.FontSize = 14;
    % grid on
    % yscale log
    % ax.YTick = 10.^(0:11); %[10^2 10^3 10^4 10^5 10^6 10^7 10^8];
    % xlim([total_time_steps-lower_value total_time_steps+upper_value])
    % ax = gca;
    % ax.FontSize = 14;
    % ax.XTick = total_time_steps-lower_value+1:20:total_time_steps+upper_value;
    % ax.XTickLabel = -lower_value:20:upper_value;
    % % ylim([min(sum(mean_x_count(total_time_steps-lower_value:total_time_steps+upper_value,:),2)) max(sum(mean_x_count(total_time_steps-lower_value:total_time_steps+upper_value,:),2))])
    % % ylim([min(mean_x_count_analy,[],"all") max(mean_x_count_analy,[],"all")])
    % ylim([10^4 0.5*10^8])
    % xlabel('time')

%% Effect of differences in adapttaion cost of the two clones on advantage of quick adaptation with smaller adaptation cost
intermediate_env = 0;
% analysis_time = [24 24*14];
analysis_time = [24 24*7];
N = 40;
total_time_steps = sum(analysis_time);

high_cost_set = [0.0495/2 0.0495*(3/4) 0.0495];
delta_cost = 0.5;%0:0.05:1; % upto 90% reduction in cost
% high_cost_set = [0.0495];
% delta_cost = 0.5; % upto 90% reduction in cost

delta_pop_per = 0.02:0.02:5;
end_env_pair = [0 1];
benefit_time_series = zeros(analysis_time(2),length(delta_cost),length(high_cost_set),length(delta_pop_per),2);
mean_cell_count = zeros(total_time_steps,length(delta_cost),length(high_cost_set),length(delta_pop_per),2);
for high_cost_indx = 1:length(high_cost_set)
    for delta_pop_per_indx = 1:length(delta_pop_per)
        for delta_indx = 1:length(delta_cost)
            for sce_indx = 1:2
                for N_indx = 1:length(N)
                    switch sce_indx
                        case 1
                            % c0_set = [0 0.0248    0.0495];
                            c0_set = [high_cost_set(high_cost_indx)*(1-delta_cost(delta_indx))    high_cost_set(high_cost_indx)];
                            x = zeros(analysis_time(2),2*length(c0_set)+1);
                            % x([1 2],:) = [repmat([0.01 0 0.04 0 0.95 0],2,1) zeros(2,1)];
                            x([1 2],:) = [repmat([delta_pop_per(delta_pop_per_indx)/100 0 1-(delta_pop_per(delta_pop_per_indx)/100) 0],2,1) zeros(2,1)];

                        case 2
                            % c0_set = [0 0.0248];
                            c0_set = [high_cost_set(high_cost_indx)*(1-delta_cost(delta_indx))];

                            x = zeros(analysis_time(2),2*length(c0_set)+1);
                            % x([1 2],:) = [repmat([0.2 0 0.8 0],2,1) zeros(2,1)];
                            x([1 2],:) = [repmat([1 0],2,1) zeros(2,1)];

                    end

                    primining_duration = 0;
                    prob_state_1_start_0 = zeros(analysis_time(2)+1,length(c0_set));
                    P_s = zeros(2,2,analysis_time(2)+2,length(c0_set));
                    for c0_indx = 1:length(c0_set)
                        q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
                        q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
                        % analytical calculation of change in cell state and benefit for two step
                        % change in the environment
                        [prob_state_1_start_0(:,c0_indx),prob_Sn_Sh_Sn_1_Sh_start_0,prob_Sn_Sl_Sn_1_Sl_start_0] = two_step_env_switch_cell_state_dyn(N(N_indx),q01,q10,end_env_pair,intermediate_env,analysis_time(2),primining_duration);
                        P_s(1,1,:,c0_indx) = prob_Sn_Sl_Sn_1_Sl_start_0;
                        P_s(1,2,:,c0_indx) = 1-P_s(1,1,:,c0_indx);
                        P_s(2,2,:,c0_indx) = prob_Sn_Sh_Sn_1_Sh_start_0;
                        P_s(2,1,:,c0_indx) = 1-P_s(2,2,:,c0_indx);
                    end

                    % Dynamics of population cell fraction and cellular benefit
                    fd = zeros(analysis_time(2),1);

                    % x([1 2],:) = [repmat([1/length(c0_set) 0],2,length(c0_set)) zeros(2,1)];
                    % x([1 2],:) = [repmat([0.01 0 0.95 0],2,1) zeros(2,1)];
                    % x([1 2],:) = [repmat([0.01 0 0.04 0 0.95 0],2,1) zeros(2,1)];
                    % x([1 2],:) = [repmat([0.2 0 0.8 0],2,1) zeros(2,1)];

                    benefit_time_series(1,delta_indx,high_cost_indx,delta_pop_per_indx,sce_indx) = reward(2);
                    for n = 2:analysis_time(2)-1
                        p = end_env_pair(2);
                        transiton_matrix = zeros(2*length(c0_set)+1);
                        growth_matrix = zeros(2*length(c0_set)+1);

                        S_l_fit = zeros(1,length(c0_set));
                        S_h_fit = zeros(1,length(c0_set));
                        temp_benefit = 0;
                        for c0_indx = 1:length(c0_set)

                            if((1-prob_state_1_start_0(n,c0_indx))~=0)
                                S_l_fit(c0_indx) = reward(2)*(1-p) +cost(2)*(p)-prob_state_1_start_0(n-1,c0_indx)*P_s(2,1,n,c0_indx)*c0_set(c0_indx)/(1-prob_state_1_start_0(n,c0_indx));
                            else
                                S_l_fit(c0_indx) = reward(2)*(1-p) +cost(2)*(p);
                            end

                            if((prob_state_1_start_0(n,c0_indx))~=0)
                                S_h_fit(c0_indx) = reward(1)*(p)+cost(1)*(1-p)-(1-prob_state_1_start_0(n-1,c0_indx))*P_s(1,2,n,c0_indx)*c0_set(c0_indx)/prob_state_1_start_0(n,c0_indx);
                            else
                                S_h_fit(c0_indx) = reward(1)*(p)+cost(1)*(1-p);
                            end

                            temp_benefit = temp_benefit + (S_l_fit(c0_indx)*(x(n,2*(c0_indx)-1)/sum(x(n,1:end-1))) + S_h_fit(c0_indx)*(x(n,2*(c0_indx))/sum(x(n,1:end-1))));

                            transiton_matrix(2*(c0_indx)-1:2*(c0_indx), 2*(length(c0_set))+1) = [-S_l_fit(c0_indx)*(1-double(S_l_fit(c0_indx)>0)); -S_h_fit(c0_indx)*(1-double(S_h_fit(c0_indx)>0))];
                            transiton_matrix(2*(c0_indx)-1:2*(c0_indx),2*(c0_indx)-1:2*(c0_indx)) = (1-transiton_matrix(2*(c0_indx)-1:2*(c0_indx), 2*(length(c0_set))+1)).*[P_s(1,1,n,c0_indx) P_s(1,2,n,c0_indx); P_s(2,1,n,c0_indx) P_s(2,2,n,c0_indx)];
                            growth_matrix(2*(c0_indx)-1:2*(c0_indx),2*(c0_indx)-1:2*(c0_indx)) = [S_l_fit(c0_indx)*double(S_l_fit(c0_indx)>0) 0; 0 S_h_fit(c0_indx)*double(S_h_fit(c0_indx)>0)];
                        end
                        benefit_time_series(n,delta_indx,high_cost_indx,delta_pop_per_indx,sce_indx) = temp_benefit;
                        transiton_matrix(end,:) = [zeros(1,2*length(c0_set)) 1];

                        x(n+1,:) = x(n,:)*(transiton_matrix+growth_matrix);
                        fd(n) = x(n,:)*sum(transiton_matrix+growth_matrix,2);
                        x(n+1,:) = x(n+1,:)/fd(n);
                    end

                    %%%%% population dynamics %%%%%%%%%
                    x = zeros(sum(analysis_time),1);
                    x(1) = 10^4;

                    for n = 1:sum(analysis_time)-1
                        if(n<=analysis_time(1))
                            f = reward(2);
                        else
                            f = (benefit_time_series(n-analysis_time(1)+1,delta_indx,high_cost_indx,delta_pop_per_indx,sce_indx));
                        end
                        if(f >= 0)
                            num_div = f*x(n);
                            num_death = 0;
                        else
                            num_div = 0;
                            num_death = -f*x(n);
                        end

                        x(n+1) = x(n) + num_div - num_death;
                        if(x(n+1) < 0)
                            x(n+1) = 0;
                        end
                    end

                    mean_cell_count(:,delta_indx,high_cost_indx,delta_pop_per_indx,sce_indx) = x;
                    % x_frac = x./sum(x,2);
                    % mean_cell_state_frac = x_frac.*[zeros(total_time_steps,length(c0_set)); cell_state_time_series(1:total_time_steps,:,prime_indx,intermediate_env_indx)];
                    % mean_x_frac = x_frac;
                end
            end
        end
    end
end
% figure('Position',[680 535 390 343])
% for N_indx = 1%:length(N)
%     hold on
%     % plot(N,mean_cell_count(end,:,2)./mean_cell_count(end,:,1),'LineWidth',2);
%     plot(1:total_time_steps,mean_cell_count(:,N_indx,2),'LineWidth',2);
%     % hold on
%     plot(1:total_time_steps,mean_cell_count(:,N_indx,1),'LineWidth',2);
% end
% yscale log
% xlim([130 300])
% figure('Position',[680 535 390 343])
% plot(N,mean_cell_count(end,:,2)./mean_cell_count(end,:,1),'LineWidth',2);
% ax = gca;
% ax.FontSize = 14;
% grid on
% xlabel('memory size (m)');
% ylabel([{'Ratio total cell counts '}, {'- enriched to mixed pop'}]);
% xlim([min(N) max(N)])

% figure('Position',[211 356 1707 372])
% fd_all = log10(mean_cell_count(end,:,:,:,2)./mean_cell_count(end,:,:,:,1));
% for high_cost_indx = 1:length(high_cost_set)
% subplot(1,3,high_cost_indx)
% fd = log10(mean_cell_count(end,:,high_cost_indx,:,2)./mean_cell_count(end,:,high_cost_indx,:,1));
% heatmap(delta_pop_per,delta_cost*100,reshape(fd,length(delta_cost),length(delta_pop_per)),"ColorLimits",[min(fd_all,[],'all') max(fd_all,[],'all')],"GridVisible","off","ColorScaling","log")
% colormap parula
% ax = gca;
% ax.FontSize = 12;
% ylabel('% cost difference');
% xlabel('low cost cell fraction');
% end

% color_order = ["#0072BD", "#D95319","#EDB120","#77AC30"];
% figure('Position',[1069 555 466 251])
% for high_cost_indx = 1:length(high_cost_set)
% plot(delta_cost*100,mean_cell_count(end,:,high_cost_indx,1,1),'-','LineWidth',2,Color=color_order(high_cost_indx));
% hold on
% plot(delta_cost*100,mean_cell_count(end,:,high_cost_indx,1,2),'-.','LineWidth',2,Color=color_order(high_cost_indx));
% end
% yscale log
% ylabel([{'End time point'},{'total cell counts'}]);
% ax = gca;
% ax.YTick = 10.^(2:8);
% ax.FontSize = 14;
% grid on
% xlabel('% cost difference');
% % xlim([min(N) max(N)])
% 
% figure('Position',[1069 555 466 251])
% for high_cost_indx = 1:length(high_cost_set)
% plot(delta_cost*100,log10(mean_cell_count(end,:,high_cost_indx,1,2)./mean_cell_count(end,:,high_cost_indx,1,1)),'LineWidth',2);
% hold on
% end
% % yscale log
% ylabel(['log_{10}(fold change)']);
% ax = gca;
% % ax.YTick = 0:25:200;
% ax.FontSize = 14;
% grid on
% xlabel('% cost difference');
% % xlim([min(N) max(N)])

color_order = ["#0072BD", "#D95319","#EDB120","#77AC30"];
figure('Position',[1069 555 466 251])
for high_cost_indx = 1:length(high_cost_set)
plot(delta_pop_per,reshape(mean_cell_count(end,1,high_cost_indx,:,1),1,[]),'-','LineWidth',2,Color=color_order(high_cost_indx));
hold on
% plot(delta_pop_per,reshape(mean_cell_count(end,1,high_cost_indx,:,2),1,[]),'-.','LineWidth',2,Color=color_order(high_cost_indx));
end
yscale log
ylabel([{'End time point'},{'total cell counts'}]);
ax = gca;
ax.YTick = 10.^(2:8);
ax.FontSize = 14;
grid on
xlabel('percentage of low-cost subpopulation');
% xlim([min(N) max(N)])

figure('Position',[1069 555 466 251])
for high_cost_indx = 1:length(high_cost_set)
plot(delta_pop_per,reshape(log10(mean_cell_count(end,1,high_cost_indx,:,2)./mean_cell_count(end,1,high_cost_indx,:,1)),1,[]),'LineWidth',2);
hold on
end
ylabel(['log_{10}(fold change)']);
ax = gca;
ax.FontSize = 14;
grid on
% yscale log
xlabel('percentage of low-cost subpopulation');
% xlim([min(N) max(N)])