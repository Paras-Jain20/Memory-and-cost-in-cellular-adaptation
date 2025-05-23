% close all

results = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\Phenotypic_adaptation_Costly_transition\Figures for manuscript';

state_label = {'l', 'h'};
model_label = {'Undated', 'Dated'};

% Periodic environment
handles = feval(@Phenotypic_adaptation_models);
total_time_steps = 5000000;

N = 9; % memory capacity
q = 0.5; % indifference environment zero adaptation cost
p = 0; % environmental state; value doesn't matter but needs to be present
num_rep= 3; % number of replicates of numerical simulation

T_set = 10; %2:2:20;%[1:10 15:5:100]; % time period

alpha = 1; % difference between (r_h^+)-(r_l^-)
beta = alpha*q./(1-q);  % indifference environment (p_I)
cost = [1 1]; % (r_h^-), (r_l^-)
reward = [alpha+cost(2) beta+cost(1)]; % (r_h^+), (r_l^+)

% c0 = 0:min(alpha,beta(q_indx))/5:(min(alpha,beta(q_indx))-min(alpha,beta(q_indx))/100);
c0 = (min(alpha,beta)-min(alpha,beta)/2);

q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

% env = binornd(1,p,[total_time_steps+N,1]); % the state of the environment along each axis and timesteps

num_models = 2;

%%
num_benefit_centered_moments = zeros(2,length(T_set),num_models,num_rep); % first two indices corresponds to the number of p and model values
num_resi_time_centered_moments = zeros(4,2,length(T_set),num_models,num_rep); % last two indices corresponds to the number of p and model values
cell_state = zeros(total_time_steps,length(T_set),num_models,num_rep);
pi_est = zeros(total_time_steps,length(T_set),num_models,num_rep);
num_joint_resi_time_centered_moments = zeros(6,length(T_set),num_models,num_rep);


% num_mean_resi_frac= zeros(4,2); % last two indices corresponds to the number of p and model values
T_indx = 1;
for T = T_set

    env = discrete_period_env(T,total_time_steps+N); % the state of the environment along each axis and timesteps

    for model_indx = 2%1:2
        switch model_indx
            case 1
                temp_cell_state = zeros(total_time_steps,num_rep);
                temp_pi_est = zeros(total_time_steps,num_rep);
                for rep_indx = 1:num_rep
                    disp(['running analyses for model ' num2str(model_indx) ' T ' num2str(T) ' rep ' num2str(rep_indx)])

                    [temp_cell_state(:,rep_indx),benefit,avg_benefit,avg_benefit_with_time,temp_pi_est(:,rep_indx)] = handles{2}(reward,cost,total_time_steps, env, N, p, c0);
                    num_benefit_centered_moments(:,T_indx,model_indx,rep_indx) = [avg_benefit; var(benefit)];

                    [residence_time_state_0,residence_time_state_1] = residence_time_cal_cell_state_based(temp_cell_state(:,rep_indx));
                    if(length(residence_time_state_0) >= 3 || length(residence_time_state_1) >= 3)
                        num_resi_time_centered_moments(:,:,T_indx,model_indx,rep_indx) = [mean(residence_time_state_0) mean(residence_time_state_1); var(residence_time_state_0) var(residence_time_state_1); skewness(residence_time_state_0) skewness(residence_time_state_1); kurtosis(residence_time_state_0) kurtosis(residence_time_state_1)];
                        min_len = min(length(residence_time_state_0),length(residence_time_state_1));
                        num_joint_resi_time_centered_moments(:,T_indx,model_indx,rep_indx) = [mean(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) ...
                            mean(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len)))...
                            mean((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len)))];

                    else
                        num_resi_time_centered_moments(:,:,T_indx,model_indx,rep_indx) = [mean(residence_time_state_0) mean(residence_time_state_1); 0 0; 0 0; 0 0];
                        min_len = min(length(residence_time_state_0),length(residence_time_state_1));
                        num_joint_resi_time_centered_moments(:,T_indx,model_indx,rep_indx) = [mean(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) 0 ...
                            mean(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) 0 ...
                            mean((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) 0];
                    end
                end
                cell_state(:,T_indx,model_indx,:) = temp_cell_state;
                pi_est(:,T_indx,model_indx,:) = temp_pi_est;

            case 2

                temp_cell_state = zeros(total_time_steps,num_rep);
                temp_pi_est = zeros(total_time_steps,num_rep);
                for rep_indx = 1:num_rep
                    disp(['running analyses for model ' num2str(model_indx) ' T ' num2str(T) ' rep ' num2str(rep_indx)])

                    [temp_cell_state(:,rep_indx),benefit,avg_benefit,avg_benefit_with_time,temp_pi_est(:,rep_indx)] = handles{1}(reward,cost,total_time_steps, env, N, p, c0);
                    num_benefit_centered_moments(:,T_indx,model_indx,rep_indx) = [avg_benefit; var(benefit)];

                    [residence_time_state_0,residence_time_state_1] = residence_time_cal_cell_state_based(temp_cell_state(:,rep_indx));

                    if(length(residence_time_state_0) >= 3 || length(residence_time_state_1) >= 3)
                        num_resi_time_centered_moments(:,:,T_indx,model_indx,rep_indx) = [mean(residence_time_state_0) mean(residence_time_state_1); var(residence_time_state_0) var(residence_time_state_1); skewness(residence_time_state_0) skewness(residence_time_state_1); kurtosis(residence_time_state_0) kurtosis(residence_time_state_1)];
                        min_len = min(length(residence_time_state_0),length(residence_time_state_1));
                        num_joint_resi_time_centered_moments(:,T_indx,model_indx,rep_indx) = [mean(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) ...
                            mean(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len)))...
                            mean((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) var((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len)))];

                    else
                        num_resi_time_centered_moments(:,:,T_indx,model_indx,rep_indx) = [mean(residence_time_state_0) mean(residence_time_state_1); 0 0; 0 0; 0 0];
                        min_len = min(length(residence_time_state_0),length(residence_time_state_1));
                        num_joint_resi_time_centered_moments(:,T_indx,model_indx,rep_indx) = [mean(residence_time_state_0(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) 0 ...
                            mean(residence_time_state_1(1:min_len)./(residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) 0 ...
                            mean((residence_time_state_0(1:min_len)+residence_time_state_1(1:min_len))) 0];
                    end

                end
                cell_state(:,T_indx,model_indx,:) = temp_cell_state;
                pi_est(:,T_indx,model_indx,:) = temp_pi_est;
        end
    end
    T_indx = T_indx + 1;

end

figure
histogram(residence_time_state_0,Normalization="probability",BinWidth=1)
xlim([0 100])
figure
histogram(residence_time_state_1,Normalization="probability",BinWidth=1)
xlim([0 100])

env = discrete_period_env(50,total_time_steps+N); % the state of the environment along each axis and timesteps

figure
plot(1:total_time_steps,env(N+1:end))
hold on
plot(1:total_time_steps,cell_state(:,10,1,1))
xlim([0 500])

% save(['Data periodic env model comparsion N ' num2str(N) ' q ' num2str(q) ' c0 ' num2str(c0) ' p2 ' num2str(p2) '.mat'])

% num_resi_time_centered_moments = num_resi_time_centered_moments(num_resi_time_centered_moments ~= 0);
cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\Phenotypic_adaptation_Costly_transition\Figures for Work Presentation')

cd(results)

save(['Data periodic env model comparsion N ' num2str(N) ' q ' num2str(q) ' c0 ' num2str(c0) '.mat'])

%%
c0 = 0; %0.21429
load(['Data periodic env model comparsion N ' num2str(N) ' q ' num2str(q) ' c0 ' num2str(c0) '.mat']);

% Figure 1

% line style and color code for type of result (analytical/numerical) and memory models
line_style = {'-', '--'}; % for analytical vs numerical result
line_color = ["#0072BD","#D95319"];

% subplot  1 - comparing central moments of residence time in cell state 0 and 1

% color_order = ["#0072BD","#4DBEEE", "#D95319", "#DD8E6A"];

figure('Units','normalized','Position',[0 0 1 1])
for state = 1:2
    subplot(2,3,state)
    y = zeros(length(T_set),2,2);
    for mom = 1:2
        % subplot(2,4,(state-1)*4+mom)

        for model_indx = 1:2
            % numerical
            for T_indx = 1:length(T_set)
                temp_num_resi_time_centered_moments = reshape(num_resi_time_centered_moments(mom,state,T_indx,model_indx,:),1,num_rep);
                non_zeros_indices = find(temp_num_resi_time_centered_moments ~=0);
                y(T_indx,model_indx,mom) = (mean(temp_num_resi_time_centered_moments(non_zeros_indices)));
            end
        end
    end

    % errorbar(T_set,y(:,:,1),y(:,:,2));
    b = bar(categorical(T_set),y(:,:,1));

    hold on
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(y(:,:,1));
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',y(:,:,1),sqrt(y(:,:,2)),'k','linestyle','none','LineWidth',1.25);
    hold off

    set(gca,'YScale','log')

    ax = gca;
    ax.FontSize = 18;
    grid on
    xlabel('Time period')
    %         colororder(color_order)

    ylabel(['\tau_' state_label{state}])

end
legend('Undated', 'Dated')
saveas(gcf,['Comparing Residence Time Central Moments periodic env q ' num2str(q) ' N ' num2str(N) ' c0 ' num2str(c0) '.png']);

% subplot  2 - comparing residence time fractions - E(x)/(E(x) + E(y)) and
% summary statistics of x/(x+y) where x and y's are residence times

% fraction mean residence time
figure('Units','normalized','Position',[0 0 1 1])
subplot(2,3,1)

y_1_frac = zeros(length(T_set),4);
for model_indx = 1:2

    % numerical results
    num_y_1_frac = zeros(length(T_set),1);

    for T_indx = 1:length(T_set)

        temp_num_resi_time_centered_moments_state_0 = reshape(num_resi_time_centered_moments(1,1,T_indx,model_indx,:),1,num_rep);
        non_zeros_indices = find(temp_num_resi_time_centered_moments_state_0 ~=0);
        temp_num_resi_time_centered_moments_state_1 = reshape(num_resi_time_centered_moments(1,2,T_indx,model_indx,:),1,num_rep);

        if(~isempty(non_zeros_indices))
            num_y_1_frac(T_indx) = (mean(temp_num_resi_time_centered_moments_state_1(non_zeros_indices)))./sum((mean(temp_num_resi_time_centered_moments_state_0(non_zeros_indices)))+(mean(temp_num_resi_time_centered_moments_state_1(non_zeros_indices))));
        else
            num_y_1_frac(T_indx) = 0;
        end
    end
    y_1_frac(:,model_indx) = num_y_1_frac;

end
% y_1_frac(1:3,2) =ones(3,1);

bar(categorical(T_set),y_1_frac,'BarWidth',.8);
% colororder(color_order)
ax = gca;
ax.FontSize = 18;
grid on
xlabel('Time period')
ylabel('fraction mean residence time')
legend('Undated', 'Dated')

saveas(gcf,['Comparing Mean Residence time fraction periodic env q ' num2str(q) ' N ' num2str(N) ' c0 ' num2str(c0) '.png']);

% mean fraction residence times

figure('Units','normalized','Position',[0 0 1 1])
for state = 1:2
    subplot(2,3,state)

    y = zeros(length(T_set),2,2);
    for mom = 1:2

        for model_indx = 1:2

            % numerical
            for T_indx = 1:length(T_set)
                temp_num_resi_time_centered_moments = reshape(num_joint_resi_time_centered_moments((state-1)*2+1+(mom-1),T_indx,model_indx,:),1,num_rep);
                non_zeros_indices = find(temp_num_resi_time_centered_moments ~=0);
                y(T_indx,model_indx,mom) = (mean(temp_num_resi_time_centered_moments(non_zeros_indices)));
            end

        end
    end
    b = bar(categorical(T_set),y(:,:,1));
    hold on
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(y(:,:,1));
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',y(:,:,1),sqrt(y(:,:,2)),'k','linestyle','none','LineWidth',1.25);
    hold off

    %         colororder(color_order)
    ax = gca;
    ax.FontSize = 18;
    grid on
    xlabel('Time period')

    ylabel(['Fraction residence time'])
end
legend('Undated', 'Dated')
saveas(gcf,['Comparing One Cycle Residence Time fraction Central Moments period env q ' num2str(q) ' N ' num2str(N) ' c0 ' num2str(c0) '.png']);


% summary statistics of total cycle time

figure('Units','normalized','Position',[0 0 1 1])
subplot(2,3,1)
y = zeros(length(T_set),2,2);
for mom = 1:2


    for model_indx = 1:2

        % % numerical
        for T_indx = 1:length(T_set)
            temp_num_resi_time_centered_moments = reshape(num_joint_resi_time_centered_moments(4+mom,T_indx,model_indx,:),1,num_rep);
            non_zeros_indices = find(temp_num_resi_time_centered_moments ~=0);
            y(T_indx,model_indx,mom) = (mean(temp_num_resi_time_centered_moments(non_zeros_indices)));
        end

    end
end
b = bar(T_set,y(:,:,1));

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(y(:,:,1));
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',y(:,:,1),sqrt(y(:,:,2)),'k','linestyle','none','LineWidth',1.25);
hold off

set(gca,'YScale','log')
ylim([1 max(y,[],'all')+10])
%         colororder(color_order)
ax = gca;
ax.FontSize = 18;
grid on
xlabel('Time period')
ylabel(['One cycle duration'])
legend('Undated', 'Dated')
saveas(gcf,['Comparing One Cycle length Central Moments periodic env q ' num2str(q) ' N ' num2str(N) ' c0 ' num2str(c0) '.png']);

% comparing central momemnts of benefit between two models
% color_order = ["#0072BD","#4DBEEE", "#D95319", "#DD8E6A"];
figure('Units','normalized','Position',[0 0 1 1])
% for state = 1:2
subplot(2,3,1)
for mom = 1:2
    
    y(:,:,mom) = [reshape(mean(num_benefit_centered_moments(mom,:,1,:),4),length(T_set),1) reshape(mean(num_benefit_centered_moments(mom,:,2,:),4),length(T_set),1)];
    % err_y = [zeros(4,1) reshape(var(num_resi_time_centered_moments(2,1,:,1,:),0,5),length(T_set),1) zeros(4,1) reshape(var(num_resi_time_centered_moments(2,1,:,2,:),0,5),length(T_set),1)];
end

b = bar(categorical(T_set),y(:,:,1));
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(y(:,:,1));
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',y(:,:,1),sqrt(y(:,:,2)),'k','linestyle','none','LineWidth',1.25);
hold off
% errorbar(categorical(T_set),y,err_y);
%     colororder(color_order)
ax = gca;
ax.FontSize = 18;
grid on
xlabel('Time period')
ylabel(['R_n'])


legend('Undated', 'Dated')
saveas(gcf,['Comparing Cell Benefit Central Moments periodic env q ' num2str(q) ' N ' num2str(N) ' c0 ' num2str(c0) '.png']);

% color_order = ["#0072BD", "#D95319", "#EDB120","#7E2F8E"];
% for plotting autocorrelation of cell memory and state for the two models across p values
marker_type = {'+','o'};
for model_indx = 1:2
    figure('Units','normalized','Position',[0 0 1 1]);

    for T_indx = 2:2:length(T_set)
        autocorr_memory = zeros(num_rep,2*T_set(end)+1);
        autocorr_cell_state = zeros(num_rep,2*T_set(end)+1);
        for rep_indx = 1:num_rep
            autocorr_memory(rep_indx,:) = autocorr(pi_est(:,T_indx,model_indx,rep_indx),NumLags=2*T_set(end));
            autocorr_cell_state(rep_indx,:) = autocorr(cell_state(:,T_indx,model_indx,rep_indx),NumLags=2*T_set(end));
        end
        subplot(2,2,1)
        plot(0:size(autocorr_memory,2)-1, mean(autocorr_memory,1),'LineWidth',2);
        hold on
        grid on
        xlabel('Lags')
        ylabel('Autocorrelation (cell memory)')
        ax = gca;
        ax.FontSize = 14;
        ax.YTick = -1:0.2:1;
        ylim([-1 1])
        %         colororder(color_order);

        subplot(2,2,2)
        plot(0:size(autocorr_cell_state,2)-1, mean(autocorr_cell_state,1),'LineWidth',2);
        hold on
        grid on
        xlabel('Lags')
        ylabel('Autocorrelation (cell state)')
        ax = gca;
        ax.FontSize = 14;
        ax.YTick = -1:0.2:1;
        ylim([-1 1])
        %         colororder(color_order);
    end
    lg = legend(num2str(T_set(2:2:end)'));
    title(lg,'Time period')
    saveas(gcf,['Numerical Autocorrelation of memory and cell state model ' model_label{model_indx} ' periodic env q ' num2str(q) ' N ' num2str(N) ' c0 ' num2str(c0) '.png']);
end
