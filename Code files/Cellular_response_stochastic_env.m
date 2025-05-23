close all

%%
% Directories
% results = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\Phenotypic_adaptation_Costly_transition\Figures for manuscript';

handles = feval(@Phenotypic_adaptation_models);
total_time_steps = 100000;

state_label = {'l', 'h'};
model_label = {'Undated', 'Dated'};

N = 8; % memory capacity
q = 0.3; % indifference environment (p_I)
p = 0.4; % stochastic environmental state
num_rep= 20; % number of replicate for numberical simulations

alpha = 1; % difference between (r_h^+)-(r_l^-) 
beta = alpha*q./(1-q); % (r_l^+)-(r_h^-)
cost = [1 1]; % (r_h^-), (r_l^-)
reward = [alpha+cost(2) beta+cost(1)]; % (r_h^+), (r_l^+)

c0_set = linspace(0,min(alpha,beta)-min(alpha,beta)/100,3); % adaptation cost
num_models = 2;

resi_time_centered_moments = zeros(4,2,length(c0_set),num_models); % last two indices corresponds to the number of p and model values
benefit_centered_moments = zeros(2,length(c0_set),num_models); % last two indices corresponds to the number of p and model values
joint_prob_residence_time_state_moments = zeros(num_models,length(c0_set),6); % first four columns represent Expectation and Variance of x/(x+y) and y/(x+y); next two columns represent E(x+y) and Var(x+y) 

prob_residence_time_state_0 = cell(length(c0_set),num_models);
prob_residence_time_state_1 = cell(length(c0_set),num_models);

joint_prob_flag = zeros(2,length(c0_set),num_models);
%% Analytical calculations
for model_indx = 1:2
    for c0_indx = 1:length(c0_set)
        q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
        q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

        switch model_indx
            case 1
                disp(['running analytical calculations for model ' num2str(model_indx) ' c0 ' num2str(c0_set(c0_indx))])

                P = transition_matrix_P_cal(N,p);
                resi_time_centered_moments(:,:,c0_indx,model_indx) = undated_memory_marginal_resi_time_moments(N,q01,q10,P);
                [~,~,cond_prob_Sn_Sn_1] = undated_memory_cond_prob_cell_state(N,q01,q10,P);

                [benefit_centered_moments(:,c0_indx,model_indx)] = benefit_centered_moments_cal(p,reward, cost, resi_time_centered_moments(:,:,c0_indx,model_indx),cond_prob_Sn_Sn_1,c0_set(c0_indx));

            case 2
                disp(['running analytical calculations for model ' num2str(model_indx) ' c0 ' num2str(c0_set(c0_indx))])
                resi_time_centered_moments(:,:,c0_indx,model_indx) = dated_memory_marginal_resi_time_moments(N,p,q01,q10);
                [~,~,cond_prob_Sn_Sn_1] = dated_memory_cond_prob_cell_state(N,p,q01,q10);

                [benefit_centered_moments(:,c0_indx,model_indx)] = benefit_centered_moments_cal(p,reward, cost, resi_time_centered_moments(:,:,c0_indx,model_indx),cond_prob_Sn_Sn_1,c0_set(c0_indx));
        end
    end
end

%% Numerical Calculations
c0_set_num = linspace(0,min(alpha,beta)-min(alpha,beta)/100,10);

num_benefit_centered_moments = zeros(2,length(c0_set_num),num_models,num_rep); % first two indices corresponds to the number of p and model values
num_resi_time_centered_moments = zeros(4,2,length(c0_set_num),num_models,num_rep); % last two indices corresponds to the number of p and model values
num_joint_resi_time_centered_moments = zeros(6,length(c0_set_num),num_models,num_rep);

cell_state = zeros(total_time_steps,length(c0_set_num),num_models,num_rep);
pi_est = zeros(total_time_steps,length(c0_set_num),num_models,num_rep);

for model_indx = 1:2
    c0_indx = 1;
    for c0 = c0_set_num
        q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
        q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

        switch model_indx
            case 1
                temp_cell_state = zeros(total_time_steps,num_rep);
                temp_pi_est = zeros(total_time_steps,num_rep);
                for rep_indx = 1:num_rep
                    disp(['running numerical analyses for model ' num2str(model_indx) ' c0 ' num2str(c0) ' rep ' num2str(rep_indx)])
                    env = binornd(1,p,[total_time_steps+N,1]); % the state of the environment along each axis and timesteps

                    [temp_cell_state(:,rep_indx),benefit,avg_benefit,avg_benefit_with_time,temp_pi_est(:,rep_indx)] = handles{2}(reward,cost,total_time_steps, env, N, p, c0);
                    num_benefit_centered_moments(:,c0_indx,model_indx,rep_indx) = [avg_benefit; var(benefit)];

                    [residence_time_state_0,residence_time_state_1] = residence_time_cal_cell_state_based(temp_cell_state(:,rep_indx));
                    if(length(residence_time_state_0) >= 3 || length(residence_time_state_1) >= 3)
                        num_resi_time_centered_moments(:,:,c0_indx,model_indx,rep_indx) = [mean(residence_time_state_0) mean(residence_time_state_1); var(residence_time_state_0) var(residence_time_state_1); skewness(residence_time_state_0) skewness(residence_time_state_1); kurtosis(residence_time_state_0) kurtosis(residence_time_state_1)];
                        min_len = min(length(residence_time_state_0),length(residence_time_state_1));
                        num_joint_resi_time_centered_moments(:,c0_indx,model_indx,rep_indx) = [mean(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) ...
                            mean(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len)))...
                            mean((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len)))];

                    else
                        num_resi_time_centered_moments(:,:,c0_indx,model_indx,rep_indx) = zeros(4,2);
                        num_joint_resi_time_centered_moments(:,c0_indx,model_indx,rep_indx) = zeros(6,1);
                    end
                end
                cell_state(:,c0_indx,model_indx,:) = temp_cell_state;
                pi_est(:,c0_indx,model_indx,:) = temp_pi_est;

            case 2
                temp_cell_state = zeros(total_time_steps,num_rep);
                temp_pi_est = zeros(total_time_steps,num_rep);
                for rep_indx = 1:num_rep
                    disp(['running numerical analyses for model ' num2str(model_indx) ' c0 ' num2str(c0) ' rep ' num2str(rep_indx)])
                    env = binornd(1,p,[total_time_steps+N,1]); % the state of the environment along each axis and timesteps

                    [temp_cell_state(:,rep_indx),benefit,avg_benefit,avg_benefit_with_time,temp_pi_est(:,rep_indx)] = handles{1}(reward,cost,total_time_steps, env, N, p, c0);
                    num_benefit_centered_moments(:,c0_indx,model_indx,rep_indx) = [avg_benefit; var(benefit)];

                    [residence_time_state_0,residence_time_state_1] = residence_time_cal_cell_state_based(temp_cell_state(:,rep_indx));
                    if(length(residence_time_state_0) >= 3 || length(residence_time_state_1) >= 3)
                        num_resi_time_centered_moments(:,:,c0_indx,model_indx,rep_indx) = [mean(residence_time_state_0) mean(residence_time_state_1); var(residence_time_state_0) var(residence_time_state_1); skewness(residence_time_state_0) skewness(residence_time_state_1); kurtosis(residence_time_state_0) kurtosis(residence_time_state_1)];
                        min_len = min(length(residence_time_state_0),length(residence_time_state_1));
                        num_joint_resi_time_centered_moments(:,c0_indx,model_indx,rep_indx) = [mean(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) ...
                            mean(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len)))...
                            mean((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len)))];

                    else
                        num_resi_time_centered_moments(:,:,c0_indx,model_indx,rep_indx) = zeros(4,2);
                        num_joint_resi_time_centered_moments(:,c0_indx,model_indx,rep_indx) = zeros(6,1);
                    end
                end
                cell_state(:,c0_indx,model_indx,:) = temp_cell_state;
                pi_est(:,c0_indx,model_indx,:) = temp_pi_est;
        end
        c0_indx = c0_indx + 1;
    end
end

% Saving workspace
% cd(results)
% save(['Data stochastic environment model comparsion N ' num2str(N) ' q ' num2str(q) ' p ' num2str(p) '.mat'])

%% Loading workspace
% load(['Data stochastic environment model comparsion N ' num2str(N) ' q ' num2str(q) ' p ' num2str(p) '.mat']);

%% comparing central moments of residence time in cell state 0 and 1

color_order = ["#0072BD", "#D95319"];

figure('Units','normalized','Position',[0 0 1 1])
line_type = {'-','--'};
Marker_type = {'*','o'};
for state = 1:2
    y = zeros(length(c0_set),2,2);
    y_num = zeros(length(c0_set_num),2,2);
    subplot(2,3,1)
    for model_indx = 1%:2
        for mom = 1:2
            % analytical
            y(:,model_indx,mom) = reshape(resi_time_centered_moments(mom,state,:,model_indx),length(c0_set),1);

            % numerical
            for c0_indx = 1:length(c0_set_num)
                temp_num_resi_time_centered_moments = reshape(num_resi_time_centered_moments(mom,state,c0_indx,model_indx,:),1,num_rep);
                non_zeros_indices = find(temp_num_resi_time_centered_moments ~=0);
                y_num(c0_indx,model_indx,mom) = (mean(temp_num_resi_time_centered_moments(non_zeros_indices)));
            end
        end
        % Plot the errorbars
        plot(c0_set,y(:,model_indx,1),line_type{model_indx},'LineWidth',2,'Color',color_order(state));
        hold on
        plot(c0_set_num,y_num(:,model_indx,1),Marker_type{model_indx},'MarkerSize',12,'Color',color_order(state));
    end
end

hold on
set(gca,'YScale','log')
ax = gca;
% ax.XTickLabel = c0_set;
ax.FontSize = 16;
grid on
xlabel('Transition cost (c_{lh} = c_{hl})')
ylabel(['Residence times']);
% legend(['\tau_' state_label{1} '^A'], ['\tau_' state_label{1} '^N'],['\tau_' state_label{2} '^A'], ['\tau_' state_label{2} '^N'])
% colororder(color_order)
xlim([min(c0_set) max(c0_set)])

% legend('Undated_A','Undated_N','Dated_A','Dated_N')
saveas(gcf,['Comparing Residence Time Central Moments q ' num2str(q) ' N ' num2str(N) ' p ' num2str(p) '.png']);

%% comparing mean fraction residence times and central momemnts of benefit between two models

% fraction mean residence time
figure('Position',[680 458 467 420])
color_order = ["#000000","#ADA1A1"];
% figure
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile
y_1_frac = zeros(length(c0_set),2);
y_1_frac_num = zeros(length(c0_set_num),2);

Marker_type = {'*','o'};

for model_indx = 1%:2
    % analytical results 
    y_1_frac(:,model_indx)  = reshape(resi_time_centered_moments(1,2,:,model_indx),length(c0_set),1)./sum(reshape(resi_time_centered_moments(1,1,:,model_indx),length(c0_set),1)+reshape(resi_time_centered_moments(1,2,:,model_indx),length(c0_set),1),2);

    % numerical results

    for c0_indx = 1:length(c0_set_num)

        temp_num_resi_time_centered_moments_state_0 = reshape(num_resi_time_centered_moments(1,1,c0_indx,model_indx,:),1,num_rep);
        non_zeros_indices = find(temp_num_resi_time_centered_moments_state_0 ~=0);
        temp_num_resi_time_centered_moments_state_1 = reshape(num_resi_time_centered_moments(1,2,c0_indx,model_indx,:),1,num_rep);

        if(~isempty(non_zeros_indices))
            y_1_frac_num(c0_indx,model_indx) = (mean(temp_num_resi_time_centered_moments_state_1(non_zeros_indices)))./sum((mean(temp_num_resi_time_centered_moments_state_0(non_zeros_indices)))+(mean(temp_num_resi_time_centered_moments_state_1(non_zeros_indices))));
        else
            y_1_frac_num(c0_indx,model_indx) = 0;
        end
    end
    % yyaxis left
    plot(c0_set,y_1_frac(:,model_indx),line_type{model_indx},'LineWidth',2,'Color',color_order(model_indx));
    hold on
    plot(c0_set_num,y_1_frac_num(:,model_indx),Marker_type{model_indx},'MarkerSize',12,'Color',color_order(model_indx));
end
colororder(color_order)
ax = gca;
ax.FontSize = 16;
grid on
% xlabel('Transition cost (c_0)')
ylabel('Frac Resi Time')
xlim([min(c0_set) max(c0_set)])
ax.YColor = 'k';
% Remove only the x-axis
set(gca, 'XColor', 'none');
% legend('Undated_A','Undated_N','Dated_A','Dated_N')

%%%%%%%%%%%%%%%%%%%%

% color_order = ["#0072BD","#4DBEEE", "#D95319", "#DD8E6A"];
% figure('Units','normalized','Position',[0 0 1 1])
% % for state = 1:2
nexttile
for mom = 1:2
    y(:,:,mom) = [reshape(benefit_centered_moments(mom,:,1),length(c0_set),1) reshape(benefit_centered_moments(mom,:,2),length(c0_set),1)];
end

for mom = 1:2
    y_num(:,:,mom) = [reshape(mean(num_benefit_centered_moments(mom,:,1,:),4),length(c0_set_num),1) reshape(mean(num_benefit_centered_moments(mom,:,2,:),4),length(c0_set_num),1)];
end

% subplot(2,3,1)
% Plot the errorbars
for model_indx= 1%:2
    % yyaxis right
    plot(c0_set,y(:,model_indx,1),line_type{model_indx},'LineWidth',2,'Color',color_order(model_indx));
    hold  on
    plot(c0_set_num,y_num(:,model_indx,1),Marker_type{model_indx},'MarkerSize',12,'Color',color_order(model_indx));
end
hold off
colororder(color_order)
ax = gca;
ax.FontSize = 16;
ax.YColor = 'k';
grid on
xlabel('Transition cost (c_{lh} = c_{hl})')
ylabel('Benefit')
xlim([min(c0_set) max(c0_set)])

% legend('Undated_A','Undated_N','Dated_A','Dated_N')
saveas(gcf,['Comparing Mean residence time and Cell Benefit Central Moments q ' num2str(q) ' N ' num2str(N) ' p ' num2str(p) '.png']);

%% summary statistics of x/(x+y) where x and y's are residence times

figure('Units','normalized','Position',[0 0 1 1])
for state = 1:2
    y = zeros(length(c0_set),4,2);
    for model_indx = 1:2

    for mom = 1:2
            % analytical
            y(:,(model_indx-1)*2+1,mom) = reshape(joint_prob_residence_time_state_moments(model_indx,:,(state-1)*2+1+(mom-1)),length(c0_set),1);

            % numerical
            y_temp = zeros(length(c0_set),1);
            for c0_indx = 1:length(c0_set)
                temp_num_resi_time_centered_moments = reshape(num_joint_resi_time_centered_moments((state-1)*2+1+(mom-1),c0_indx,model_indx,:),1,num_rep);
                non_zeros_indices = find(temp_num_resi_time_centered_moments ~=0);
                y_temp(c0_indx) = (mean(temp_num_resi_time_centered_moments(non_zeros_indices)));
            end
            y(:,(model_indx-1)*2+2,mom) = y_temp;

    end
    subplot(2,3,state)
    plot(c0_set,y(:,(model_indx-1)*2+1,1),'LineWidth',2);
    hold  on
    plot(c0_set,y(:,(model_indx-1)*2+2,1),'--','LineWidth',2);
    end

    hold off
    colororder(color_order)
    ax = gca;
    ax.FontSize = 16;
    grid on
    xlabel('Transition cost (c_0)')
    ylabel(['mean fraction residence time'])
end
legend('Undated_A','Undated_N','Dated_A','Dated_N')
saveas(gcf,['Comparing One Cycle Residence Time fraction Central Moments q ' num2str(q) ' N ' num2str(N) ' p ' num2str(p) '.png']);



% %% summary statistics of total cycle time 
% 
% figure('Units','normalized','Position',[0 0 1 1])
% 
%     y = zeros(length(c0_set),4,2);
%     for mom = 1:2
%         for model_indx = 1:2
%             % analytical
%             y(:,(model_indx-1)*2+1,mom) = reshape(joint_prob_residence_time_state_moments(model_indx,:,4+mom),length(c0_set),1);
%             %             plot(categorical(c0_set),log10(y),'LineStyle',line_style{2},'Color',line_color(model_indx))
%             %             hold on;
%             % % numerical
%             y_temp = zeros(length(c0_set),1);
%             for c0_indx = 1:length(c0_set)
%                 temp_num_resi_time_centered_moments = reshape(num_joint_resi_time_centered_moments(4+mom,c0_indx,model_indx,:),1,num_rep);
%                 non_zeros_indices = find(temp_num_resi_time_centered_moments ~=0);
%                 y_temp(c0_indx) = (mean(temp_num_resi_time_centered_moments(non_zeros_indices)));
%             end
%             y(:,(model_indx-1)*2+2,mom) = y_temp;
%             % plot(categorical(c0_set),y,'LineStyle',line_style{2},'Color',line_color(model_indx));
%             % hold on;
%         end
%     end
%     subplot(2,3,1)
%     b = bar(categorical(c0_set),y(:,:,1));
%     hold on
%     % Calculate the number of groups and number of bars in each group
%     [ngroups,nbars] = size(y(:,:,1));
%     % Get the x coordinate of the bars
%     x = nan(nbars, ngroups);
%     for i = 1:nbars
%         x(i,:) = b(i).XEndPoints;
%     end
%     % Plot the errorbars
%     errorbar(x',y(:,:,1),sqrt(y(:,:,2)),'k','linestyle','none','LineWidth',1.25);
%     hold off
%     set(gca,'YScale','log')
%     ylim([1 max(y(:,:,1),[],'all')+10])
%     ylabel('One cycle duration')
%     colororder(color_order)
%     ax = gca;
%     ax.FontSize = 16;
%     grid on
%     xlabel('Env p')
%     ax.YTick = [1,10,100,1000,10000];
%     ylim([1 max(y(:,:,1),[],'all')+max(sqrt(y(:,:,2)),[],'all')]);
% legend('Undated_A','Undated_N','Dated_A','Dated_N')
% saveas(gcf,['Comparing One Cycle length Central Moments q ' num2str(q) ' N ' num2str(N) ' c0 ' num2str(c0) '.png']);
% 
% 
% %% residence time distributions 
% for model_indx = 1:2
%     figure('Units','normalized','Position',[0 0 1 1])
%     for state = 1:2
%         for c0_indx = 1:length(c0_set)
%             subplot(4,2,(c0_indx-1)*2+state)
%             if(state == 1)
%                 bar(prob_residence_time_state_0{c0_indx, model_indx});
%                 ylim([0 max([prob_residence_time_state_0{c0_indx, 1} prob_residence_time_state_0{c0_indx, 2}])])
%             else
%                 bar(prob_residence_time_state_1{c0_indx, model_indx});
%                 ylim([0 max([prob_residence_time_state_1{c0_indx, 1} prob_residence_time_state_1{c0_indx, 2}])])
%             end
%             grid on
%             title(['State S^' state_label{state} ' p ' num2str(c0_set(c0_indx))])
%             xlim([0 100])
%             ylabel('Probability')
%             xlabel(['\tau_' state_label{state} ])
%             ax =gca;
%             ax.FontSize = 12;
%         end
%     end
%     saveas(gcf,['Residence time dist model ' model_label{model_indx} ' q ' num2str(q) ' N ' num2str(N) ' c0 ' num2str(c0) '.png']);
% end
% 
% %% Autocorrelation function
% color_order = ["#0072BD", "#D95319", "#EDB120","#7E2F8E"];
% fig1 = figure('Units','normalized','Position',[0 0 1 1]);
% % for plotting autocorrelation of cell memory and state for the two models across p values
% marker_type = {'+','o'};
% for model_indx = 1:2
%     for c0_indx = 1:length(c0_set)
%         autocorr_memory = zeros(num_rep,2*N+1);
%         autocorr_cell_state = zeros(num_rep,2*N+1);
%         for rep_indx = 1:num_rep
%             autocorr_memory(rep_indx,:) = autocorr(pi_est(:,c0_indx,model_indx,rep_indx),NumLags=2*N);
%             autocorr_cell_state(rep_indx,:) = autocorr(cell_state(:,c0_indx,model_indx,rep_indx),NumLags=2*N);
%         end
%         subplot(2,2,1)
%         plot(0:size(autocorr_memory,2)-1, mean(autocorr_memory,1),'Marker',marker_type{model_indx},'MarkerSize',10,'LineWidth',2);
%         hold on
%         grid on
%         xlabel('Lags')
%         ylabel('Autocorrelation (cell memory)')
%         ax = gca;
%         ax.FontSize = 14;
%         ax.YTick = 0:0.2:1;
%         ylim([0 1])
%         colororder(color_order);
% 
% 
%         subplot(2,2,2)
%         plot(0:size(autocorr_cell_state,2)-1, mean(autocorr_cell_state,1),'Marker',marker_type{model_indx},'MarkerSize',10,'LineWidth',2);
%         hold on
%         grid on
%         xlabel('Lags')
%         ylabel('Autocorrelation (cell state)')
%         ax = gca;
%         ax.FontSize = 14;
%         ax.YTick = 0:0.2:1;
%         ylim([0 1])
%         colororder(color_order);
% 
%     end
% end
% saveas(gcf,['Numerical Autocorrelation of memory and cell state q ' num2str(q) ' N ' num2str(N) ' c0 ' num2str(c0) '.png']);