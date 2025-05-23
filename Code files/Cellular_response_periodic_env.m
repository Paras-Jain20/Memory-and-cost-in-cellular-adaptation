
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                This file contains all the codes to             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate results of cellular adaptation periodic environments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Each section headline represents the analysis carried out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
% N = 5:5:40;
% T = 5:5:40; % time period of the periodic signal
% results = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\Phenotypic_adaptation_Costly_transition\Figures for manuscript';
T = 14; % time period
N = 8; % memory size/capacity
q = 0.5; % indifference environment zero adaptation cost
p = 0; % environmental state; value doesn't matter but needs to be present
handles = feval(@Phenotypic_adaptation_models);

total_time_steps = 1000;

alpha = 1; % difference between (r_h^+)-(r_l^-) 
beta = alpha*q./(1-q); % (r_l^+)-(r_h^-)
cost = [1 1]; % (r_h^-), (r_l^-)
reward = [alpha+cost(2) beta+cost(1)]; % (r_h^+), (r_l^+)

c0 = 0;%(min(alpha,beta)-min(alpha,beta)/2); % adaptation cost

q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

%% Comparing osciallations in the estimates of env between dated and undated memory model (periodic environment)
close all
m = 2; % # of time periods to be plotted
c0 = 0;
total_time_steps = 1000;
env = discrete_period_env(T,total_time_steps+N); % the state of the environment along each axis and timesteps

% analysis for the undated memory model
P0 = transition_matrix_P_cal(N,0);
P1 = transition_matrix_P_cal(N,1);

% stationary distribution calculation
prob_analy = zeros(T,N+1);
prob_analy(end,:) = null((eye(N+1) - P0^(floor(T/2))*P1^(floor(T/2)+1))');
prob_analy(end,:)  = prob_analy(end,:)/sum(prob_analy(end,:));
for indx = 1:T-1
    if(indx <= floor(T/2))
        prob_analy(indx,:) = prob_analy(end,:)*P0^(indx);
    else
        prob_analy(indx,:) = prob_analy(end,:)*P0^(floor(T/2))*P1^(indx-floor(T/2));
    end
end


% simulation of the dated memory model
[cell_state_dated,benefit,avg_benefit,avg_benefit_with_time,pi_est_dated] = handles{1}(reward,cost,total_time_steps, env, N, p, c0);

fig1 = figure('Position',[115 407 1026 471]);
tiledlayout(2,2,'Padding','tight','TileSpacing','tight');
% plotting env
nexttile
bar(0:m*T-1,env(floor(T/2)+1:floor(T/2)+m*T),0.9,'LineWidth',1.5,'FaceColor',"#77AC30");

grid on
ax = gca;
xlim([-1 m*T])
ax.FontSize = 16;
ylabel('Env')
ax.XTick = 0:2:m*T;
ax.YTick = [0 1];
xlabel('time')

nexttile
for itr = 1:3
% simulation of the undated memory model
[cell_state_undated,benefit,avg_benefit,avg_benefit_with_time,pi_est_undated] = handles{2}(reward,cost,total_time_steps, env, N, p, c0);

% plotting pi estimate
if(N <= floor(T/2))
    hold on
    plot(0:m*T-1,pi_est_dated(1+floor(T/2)-N+1:m*T+floor(T/2)-N+1),'-','Color',"#ADA1A1",'LineWidth',2)
    plot(0:m*T-1,pi_est_undated(1+floor(T/2)-N+1:m*T+floor(T/2)-N+1),'--k','LineWidth',2)
else
    n = 2;
    while(N > n*floor(T/2))
        n = n+1;
    end
    if(mod(n,2) == 0)
        hold on
        plot(0:m*T-1,pi_est_dated(1+(n+1)*floor(T/2)-N+1:m*T+(n+1)*floor(T/2)-N+1),'-','Color',"#ADA1A1",'LineWidth',2)
        plot(0:m*T-1,pi_est_undated(1+(n+1)*floor(T/2)-N+1:m*T+(n+1)*floor(T/2)-N+1),'--k','LineWidth',2)

    else
        hold on
        plot(0:m*T-1,pi_est_dated(1+(n)*floor(T/2)-N+1:m*T+(n)*floor(T/2)-N+1),'Color','-',"#ADA1A1",'LineWidth',2)
        plot(0:m*T-1,pi_est_undated(1+(n)*floor(T/2)-N+1:m*T+(n)*floor(T/2)-N+1),'--k','LineWidth',2)
    end
end
end
plot(0:m*T,q*ones(m*T+1,1),'--r','LineWidth',1);

grid on
ax = gca;
ax.FontSize = 16;
xlabel('time')
ylabel('\pi_n')
ylim([0 1])
ax.YTick = 0:0.25:1;
ax.XTick = 0:2:m*T;
xlim([-1 m*T])

nexttile
heatmap(0:m*T-1,(N:-1:0)/N,repmat(flip([prob_analy(floor(T/2)+1:end,:); prob_analy(1:floor(T/2),:)]',1),1,2));

xlabel('time')
ylabel('\pi_n')
ax = gca;
ax.FontSize = 16;
m = 2;
pi_n_str = string((N:-1:0)/N);
for i = 1:m:N
    pi_n_str(i+1) = " ";
end
ax.YDisplayLabels= pi_n_str;

n_str = string(0:m*T-1);
for i = 1:m:m*T-1
    n_str(i+1) = " ";
end
ax.XDisplayLabels= n_str;

%     %3D bar plots for the evolution of the distribution of pi_est with
%     %time
%
%     % Create a 3D bar chart
%     nexttile
% %     figure
%     h = bar3([prob_analy(floor(T/2)+1:end,:); prob_analy(1:floor(T/2),:)]);
%
%     % % Customize the color of the bars
%     % for i = 1:length(h)
%     %     zdata = get(h(i), 'ZData');
%     %     set(h(i), 'CData', zdata, 'FaceColor', 'interp');
%     % end
%
%     % Label the axes
%     ax = gca;
%     ax.XTick = linspace(0,N,5);
%     ax.XTickLabel = linspace(0,N,5)/N;
%     ylabel('Time');
%     xlabel('Env Estimate (\pi_n)');
%     zlabel('Probability');
%     title('Dynamics of env estimate undated memory');
%     ax.FontSize = 16;

%% Using analytical calculations to compare dynamics of cell state for specific cases
close all
w = 1;
h = 1;
x_high = [];
y_high = [];
p =0;
q = 0.5;
N  = 9;
T = [3 10 30];
handles = feval(@Phenotypic_adaptation_models);

alpha = 1;
beta = alpha*q./(1-q);
cost = [1 1];
reward = [alpha+cost(2) beta+cost(1)];

c0 = 0;%(min(alpha,beta)-min(alpha,beta)/2);

q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

cum_adapt_prob =cell(length(N),length(T));
adapt_time =zeros(length(N),length(T),2);

figure('Position',[680 338 479 540])
tiledlayout(3,1,"TileSpacing","tight","Padding","tight")

for N_indx = 1:length(N)

    for T_indx = 1:length(T)

        % simulation of the dated memory model
        total_time_steps = 100;
        env = discrete_period_env(T(T_indx),total_time_steps+N(N_indx)); % the state of the environment along each axis and timesteps
        [cell_state_dated,~,avg_benefit,avg_benefit_with_time,pi_est_dated] = handles{1}(reward,cost,total_time_steps, env, N(N_indx), p, c0);

        nexttile
        m = 1;
        if(N(N_indx) <= floor(T(T_indx)/2))
            % plot(0:m*T(T_indx)-1,cell_state_dated(1+floor(T(T_indx)/2)-N(N_indx):m*T(T_indx)+floor(T(T_indx)/2)-N(N_indx)),'-','Color',"#ADA1A1",'LineWidth',2)
            bar(0:m*T(T_indx)-1,cell_state_dated(1+floor(T(T_indx)/2)-N(N_indx):m*T(T_indx)+floor(T(T_indx)/2)-N(N_indx)),'LineWidth',1.5,'FaceColor',"#ADA1A1")

        else
            n = 2;
            while(N(N_indx) > n*floor(T(T_indx)/2))
                n = n+1;
            end
            if(mod(n,2) == 0)
                  bar(0:m*T(T_indx)-1,cell_state_dated(1+(n+1)*floor(T(T_indx)/2)-N(N_indx):m*T(T_indx)+(n+1)*floor(T(T_indx)/2)-N(N_indx)),'LineWidth',1.5,'FaceColor',"#ADA1A1")

            else
                  bar(0:m*T(T_indx)-1,cell_state_dated(1+(n)*floor(T(T_indx)/2)-N(N_indx):m*T(T_indx)+(n)*floor(T(T_indx)/2)-N(N_indx)),'LineWidth',1.5,'FaceColor',"#ADA1A1")
            end
        end
        hold on

        [prob_state_1,prob_state_1_given_1] = undated_memory_cell_state_and_cond_prob_cell_state_periodic_env(N(N_indx),T(T_indx),q01,q10);

        plot(1:ceil(T(T_indx)/2),[prob_state_1(floor(T(T_indx)/2)+1:end)],'*','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2);
        plot([0 ceil(T(T_indx)/2)+1:T(T_indx)-1],[prob_state_1(floor(T(T_indx)/2)) prob_state_1(1:floor(T(T_indx)/2)-1)],'*','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2);
        % xlabel('Time Period (T)')
        % ylabel('Cell State')
        grid on
        ax = gca;
        ax.FontSize = 20;
        if(T_indx==length(T))
            xlabel('time')
        end
        title(['N ' num2str(N(N_indx)) ' T ' num2str(T(T_indx))])
if(T_indx > 1)
                ax.XTick = [0 floor(T(T_indx)/2)-1 T(T_indx)-1];
end
plot(ones(1,11)*ceil(T(T_indx)/2), 0:0.1:1,'--r','LineWidth',2);
    end
end

% cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\Phenotypic_adaptation_Costly_transition\Figures for manuscript')
% saveas(gcf,['Prob_trans_period_env q ' num2str(q) ' cost ' num2str(c0) '.png'])

%% Comparing adaptation times between dated and undated memory models in periodic environments
% Using analytical calculations to compare the prob of adaptation
% close all
w = 1;
h = 1;
x_high = [];
y_high = [];
q = 0.3;
N  = 2:12;
T = 2:40;

alpha = 1;
beta = alpha*q./(1-q);
cost = [1 1];
reward = [alpha+cost(2) beta+cost(1)];

c0 = 0;%(min(alpha,beta)-min(alpha,beta)/2);

q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

cum_adapt_prob =zeros(length(N),length(T),2);

adapt_time_dated_mem =zeros(length(N),length(T),2);
adapt_time_undated_mem =zeros(length(N),length(T),2);
adapt_time_undated_mem_std =zeros(length(N),length(T),2);

% figure('Position',[680 113 591 765]);

for N_indx = 1:length(N)

    for T_indx = 1:length(T) 

        if(N(N_indx) < T(T_indx))
            if(N(N_indx) <= floor(T(T_indx)/2))
                kmax = N(N_indx); kmin = 0;
                % for 0 -> 1 transition
                steps_for_01_trans = floor(N(N_indx)*q01)+1;
                % for 1 -> 0 transition
                steps_for_10_trans = kmax - floor(N(N_indx)*q10);
            else
                kmax = ceil(T(T_indx)/2); kmin = N(N_indx) - floor(T(T_indx)/2);
                % for 0 -> 1 transition
                steps_for_01_trans = N(N_indx) - ceil(T(T_indx)/2)+floor(N(N_indx)*q01)+1-kmin;
                % for 1 -> 0 transition
                steps_for_10_trans = N(N_indx) - floor(T(T_indx)/2) - floor(N(N_indx)*q10) + kmax;
            end
        else
            l = rem(N(N_indx),T(T_indx));
            n = (N(N_indx)-l)/T(T_indx);
            if(l <= floor(T(T_indx)/2))
                kmax = n*ceil(T(T_indx)/2)+l; kmin = n*ceil(T(T_indx)/2);
                % for 0 -> 1 transition
                steps_for_01_trans = floor(N(N_indx)*q01)+1-kmin;
                % for 1 -> 0 transition
                steps_for_10_trans = kmax - floor(N(N_indx)*q10);
            else
                kmax = (n+1)*ceil(T(T_indx)/2); kmin = n*ceil(T(T_indx)/2)+l-floor(T(T_indx)/2);
                % for 0 -> 1 transition
                steps_for_01_trans = l - ceil(T(T_indx)/2)+floor(N(N_indx)*q01)+1-kmin;
                % for 1 -> 0 transition
                steps_for_10_trans = l - floor(T(T_indx)/2) - floor(N(N_indx)*q10) + kmax;
            end
        end


        if(N(N_indx)*q01 >= kmax ||  N(N_indx)*q10 < kmin)
            adapt_time_dated_mem(N_indx,T_indx,:) = nan;
            x_high = [x_high T_indx];
            y_high = [y_high N_indx]; % Row index - 0.5
        else
            adapt_time_dated_mem(N_indx,T_indx,:) = [steps_for_01_trans/floor(T(T_indx)/2), steps_for_10_trans/ceil(T(T_indx)/2)];
        end

        centered_moments_adapt_time = undated_memory_adaptation_time_periodic_env(N(N_indx),T(T_indx),q01,q10);
        if(centered_moments_adapt_time(1,1) > floor(T(T_indx)/2))
        adapt_time_undated_mem(N_indx,T_indx,1) = nan;
        else
            adapt_time_undated_mem(N_indx,T_indx,1) = centered_moments_adapt_time(1,1)/floor(T(T_indx)/2);
        end

        if(centered_moments_adapt_time(1,2) > ceil(T(T_indx)/2))
        adapt_time_undated_mem(N_indx,T_indx,2) = nan;
        else
%             adapt_time_undated_mem(N_indx,T_indx,2) = centered_moments_adapt_time(1,2);
            adapt_time_undated_mem(N_indx,T_indx,2) = centered_moments_adapt_time(1,2)/ceil(T(T_indx)/2);

        end

        adapt_time_undated_mem_std(N_indx,T_indx,:) = (centered_moments_adapt_time(2,:)).^(0.5);
    end
end

figure('Position',[69 422 1146 523])
tiledlayout(2,2,"TileSpacing","compact","Padding","compact")
% c_lims_01 = [min([adapt_time_undated_mem(:,:,1) adapt_time_dated_mem(:,:,1)],[],'all') max([adapt_time_undated_mem(:,:,1) adapt_time_dated_mem(:,:,1)],[],'all')];
% c_lims_10 = [min([adapt_time_undated_mem(:,:,2) adapt_time_dated_mem(:,:,2)],[],'all') max([adapt_time_undated_mem(:,:,2) adapt_time_dated_mem(:,:,2)],[],'all')];

c_lims_01 = [0 1];
c_lims_10 = [0 1];


for k = 1:4
    nexttile
%     switch k
%         case 1
%             heatmap(T,N,adapt_time_dated_mem(:,:,1),'ColorLimits',c_lims_01)
%             title('0 1 transition')
%         case 2
%             heatmap(T,N,adapt_time_dated_mem(:,:,2),'ColorLimits',c_lims_10)
%             title('1 0 transition')
%         case 3
%             heatmap(T,N,adapt_time_undated_mem(:,:,1),'ColorLimits',c_lims_01)
%         case 4
%             heatmap(T,N,adapt_time_undated_mem(:,:,2),'ColorLimits',c_lims_10,'CellLabelColor','none')
%     end

    switch k
        case 1
            heatmap(T,N,adapt_time_undated_mem(:,:,1),'ColorLimits',c_lims_01)
            title('Undated')
                ylabel('Memory Size (m)')
    colorbar off

        case 2
            heatmap(T,N,adapt_time_dated_mem(:,:,1),'ColorLimits',c_lims_10)
            title('Dated')
                colorbar off

        case 3
            heatmap(T,N,adapt_time_undated_mem(:,:,2),'ColorLimits',c_lims_01)
            xlabel('Time Period (T)')
    ylabel('Memory Size (m)')
    colorbar off

        case 4
            heatmap(T,N,adapt_time_dated_mem(:,:,2),'ColorLimits',c_lims_10,'CellLabelColor','none')
            xlabel('Time Period (T)')

    end
%     xlabel('Time Period (T)')
%     ylabel('Memory Size (N)')
    ax = gca;
    ax.FontSize = 16;
    colormap parula
    count = 1;
    for i = T
        if(mod(i,4) == 0)
            ax.XDisplayLabels{count} = string(i);
        else
            ax.XDisplayLabels{count} = " ";
        end
        count = count +1;
    end

    count = 1;
    for i = N
        if(mod(i,2) == 0)
            ax.YDisplayLabels{count} = string(i);
        else
            ax.YDisplayLabels{count} = " ";
        end
        count = count +1;
    end

end

% cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\Phenotypic_adaptation_Costly_transition\Figures for manuscript')
% saveas(gcf,['Relative_Adaptation_time_period_env q ' num2str(q) ' cost ' num2str(c0) '.png'])

%% Comparing mean benefit between dated and undated memory models in periodic environments
close all

w = 1;
h = 1;
x_high = [];
y_high = [];

x_high_1 = [];
y_high_1 = [];

q = 0.3;
N  = 2:12;
T = 2:40;

% N  = 9;
% T = 3:2:40;

alpha = 1;
beta = alpha*q./(1-q);
cost = [1 1];
reward = [alpha+cost(2) beta+cost(1)];

c0 = 0;%(min(alpha,beta)-min(alpha,beta)/2);

q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

benefit = zeros(length(N),length(T),2);

% figure %('Position',[680 113 591 765]);

for N_indx = 1:length(N)

    for T_indx = 1:length(T) 

        % determining benefit for the dated memory model
        if(N(N_indx) < T(T_indx))
            if(N(N_indx) <= floor(T(T_indx)/2))
                kmax = N(N_indx); kmin = 0;
                % for 0 -> 1 transition
                steps_for_01_trans = floor(N(N_indx)*q01)+1;
                % for 1 -> 0 transition
                steps_for_10_trans = kmax - floor(N(N_indx)*q10);
            else
                kmax = ceil(T(T_indx)/2); kmin = N(N_indx) - floor(T(T_indx)/2);
                % for 0 -> 1 transition
                steps_for_01_trans = N(N_indx) - ceil(T(T_indx)/2)+floor(N(N_indx)*q01)+1-kmin;
                % for 1 -> 0 transition
                steps_for_10_trans = N(N_indx) - floor(T(T_indx)/2) - floor(N(N_indx)*q10) + kmax;
            end
        else
            l = rem(N(N_indx),T(T_indx));
            n = (N(N_indx)-l)/T(T_indx);
            if(l <= floor(T(T_indx)/2))
                kmax = n*ceil(T(T_indx)/2)+l; kmin = n*ceil(T(T_indx)/2);
                % for 0 -> 1 transition
                steps_for_01_trans = floor(N(N_indx)*q01)+1-kmin;
                % for 1 -> 0 transition
                steps_for_10_trans = kmax - floor(N(N_indx)*q10);
            else
                kmax = (n+1)*ceil(T(T_indx)/2); kmin = n*ceil(T(T_indx)/2)+l-floor(T(T_indx)/2);
                % for 0 -> 1 transition
                steps_for_01_trans = l - ceil(T(T_indx)/2)+floor(N(N_indx)*q01)+1-kmin;
                % for 1 -> 0 transition
                steps_for_10_trans = l - floor(T(T_indx)/2) - floor(N(N_indx)*q10) + kmax;
            end
        end

        if(N(N_indx)*q01 >= kmax && N(N_indx)*q10 >= kmin)
            benefit(N_indx,T_indx,1) = cost(2)*ceil(T(T_indx)/2)+ reward(2)*floor(T(T_indx)/2);
        elseif(N(N_indx)*q10 < kmin && N(N_indx)*q01 < kmax)
            benefit(N_indx,T_indx,1) = cost(1)*floor(T(T_indx)/2)+ reward(1)*ceil(T(T_indx)/2);
        elseif(N(N_indx)*q10 < kmin && N(N_indx)*q01 >= kmax) % the asymptotic average benefit is initial cell state dependent as there won't be cell state switching 
            benefit(N_indx,T_indx,1) = (cost(1)*ceil(T(T_indx)/2)+ reward(2)*floor(T(T_indx)/2)); % + 0.5*(cost(1)*floor(T/2)+ reward(1)*ceil(T/2)); assuming the initial cell state is 0
        else
            % benefit(N_indx,T_indx,1) = cost(1)*steps_for_01_trans + reward(1)*(ceil(T(T_indx)/2) - steps_for_01_trans) + cost(2)*steps_for_10_trans + reward(2)*(floor(T(T_indx)/2)-steps_for_10_trans);
            benefit(N_indx,T_indx,1) = cost(2)*steps_for_01_trans + reward(1)*(ceil(T(T_indx)/2) - steps_for_01_trans) + cost(1)*steps_for_10_trans + reward(2)*(floor(T(T_indx)/2)-steps_for_10_trans);
        end
        benefit(N_indx,T_indx,1) = benefit(N_indx,T_indx,1)/T(T_indx);
        
        if(N(N_indx)*q01 >= kmax ||  N(N_indx)*q10 < kmin)
            x_high = [x_high T_indx];
            y_high = [y_high N_indx]; % Row index - 0.5
        end


        % determining benefit for the undated memory model
        centered_moments_adapt_time = undated_memory_adaptation_time_periodic_env(N(N_indx),T(T_indx),q01,q10);
        if(centered_moments_adapt_time(1,1) > floor(T(T_indx)/2) || centered_moments_adapt_time(1,2) > ceil(T(T_indx)/2))
            x_high_1 = [x_high_1 T_indx];
            y_high_1 = [y_high_1 N_indx]; % Row index - 0.5
        end

        [prob_state_1,prob_state_1_given_1] = undated_memory_cell_state_and_cond_prob_cell_state_periodic_env(N(N_indx),T(T_indx),q01,q10);
        benefit(N_indx,T_indx,2) = benefit_centered_moments_cal_periodic_env(T(T_indx),reward,cost, prob_state_1, prob_state_1_given_1,c0); 

    end
end


% for m = 1:2
%     if (m==1)
%     plot(T,reshape(benefit(1,:,m),length(T),1),'LineWidth',2,'Color',"#ADA1A1");
%     hold on
%     else
%     plot(T,reshape(benefit(1,:,m),length(T),1),'-k','LineWidth',2);
%     end
% end
% 
% grid on
% ax =gca;
% ax.FontSize = 14;
% ax.XTick = 3:4:40;
% xlabel('Time period (T)');
% ylabel('Average Benefit');
% saveas(figure(1),['Benefit comparison period_env N ' num2str(N) 'q ' num2str(q) ' cost ' num2str(c0) '.png'])


%%% Red white blue colormap
relative_benefit = benefit(:,:,1)./benefit(:,:,2);
cmap = color_map_func();

% t1 = tiledlayout(2,1);
for j = 1:3
    figure('Position',[589 267 531 258]);
    % figure('Position',[680 629 610 249])
    nexttile
    if(j == 1)
    imagesc(benefit(:,:,1),[min(benefit(:,:,1),[],'all') max(benefit(:,:,1),[],'all')]);
    colormap parula
    elseif(j==2)
        imagesc(benefit(:,:,2),[min(benefit(:,:,2),[],'all') max(benefit(:,:,2),[],'all')]);
        colormap parula
    else
    c_lim = max([max(relative_benefit,[],'all')-1 1-min(relative_benefit,[],"all")]);
    imagesc(relative_benefit,[1-c_lim 1+c_lim]);
    colormap(cmap)
    end
    
    c = colorbar;
    % c.Limits = [0 1];
    % axis equal tight; % Adjust axis for proper grid alignment
    % Highlight specific grids (e.g., row 3, col 4)
    hold on;
    rect_color = {'w','w','k'};
    for i = 1:length(y_high)
        if(j == 1)
            rectangle('Position', [x_high(i)-0.5, y_high(i)-0.5, w, h], 'EdgeColor', rect_color{j}, 'LineWidth', 2);
        elseif(j==2)
            rectangle('Position', [x_high_1(i)-0.5, y_high_1(i)-0.5, w, h], 'EdgeColor', rect_color{j}, 'LineWidth', 2);
        end
    end

    if(j==1)
        for i = 1:length(y_high)
            rectangle('Position', [x_high(i)-0.5, y_high(i)-0.5, w, h], 'EdgeColor', rect_color{j}, 'LineWidth', 2);
        end
    elseif(j==2)
        for i = 1:length(y_high_1)
            rectangle('Position', [x_high_1(i)-0.5, y_high_1(i)-0.5, w, h], 'EdgeColor', rect_color{j}, 'LineWidth', 2);
        end
    end

    % heatmap(c0,T,benefit(:,:,j),"Colormap",cmap,"ColorLimits",[1 2],'GridVisible','off');
    ax= gca;
    ax.FontSize = 16;
    ax.YTick = [1:2:length(N)];
    ax.YTickLabels = [N(1:2:length(N))];
    ax.XTick = [3:4:length(T)];
    ax.XTickLabels = [(T(3:4:length(T)))];
    ylabel('Memory size (m)')
    xlabel('Time Period (T)')
end
% cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\Phenotypic_adaptation_Costly_transition\Figures for manuscript')
% saveas(figure(1),['Dated memory benefit period_env ' 'q ' num2str(q) ' cost ' num2str(c0) '.png'])
% saveas(figure(2),['Undated memory benefit period_env ' 'q ' num2str(q) ' cost ' num2str(c0) '.png'])
% saveas(figure(3),['Relative benefit period_env ' 'q ' num2str(q) ' cost ' num2str(c0) '.png'])

%% time dynamics of cell state with cost variation

m = 1; % number of time period to be plotted
c0 = linspace(0,(min(alpha,beta)-min(alpha,beta)/100),5);
% color_code= ["#eaeded","#bfc9ca","#95a5a6","#717d7e","#4d5656"];
fig2 = figure('Position',[589 461 1027 415]);
t1 = tiledlayout(2,2,'Padding','tight','TileSpacing','tight');

for c0_indx = 1:length(c0)

    % % simulation of the undated memory model
    % [cell_state_undated,benefit,avg_benefit,avg_benefit_with_time,pi_est_undated] = handles{2}(reward,cost,total_time_steps, env, N, p, c0(c0_indx));

    % simulation of the dated memory model
    [cell_state_dated,benefit,avg_benefit,avg_benefit_with_time,pi_est_dated] = handles{1}(reward,cost,total_time_steps, env, N, p, c0(c0_indx));

    % plotting cell state
    nexttile(t1,1)
    if(c0_indx == 2)
    if(N <= floor(T/2))
        plot(0:m*T-1,cell_state_dated(1+floor(T/2)-N+1:m*T+floor(T/2)-N+1),'--','LineWidth',2)
    else
        n = 2;
        while(N > n*floor(T/2))
            n = n+1;
        end
        if(mod(n,2) == 0)
            plot(0:m*T-1,cell_state_dated(1+(n+1)*floor(T/2)-N+1:m*T+(n+1)*floor(T/2)-N+1),'--','LineWidth',2)
        else
            plot(0:m*T-1,cell_state_dated(1+(n)*floor(T/2)-N+1:m*T+(n)*floor(T/2)-N+1),'--','LineWidth',2)
        end
    end
    else
        if(N <= floor(T/2))
            plot(0:m*T-1,cell_state_dated(1+floor(T/2)-N+1:m*T+floor(T/2)-N+1),'-','LineWidth',2)
        else
            n = 2;
            while(N > n*floor(T/2))
                n = n+1;
            end
            if(mod(n,2) == 0)
                plot(0:m*T-1,cell_state_dated(1+(n+1)*floor(T/2)-N+1:m*T+(n+1)*floor(T/2)-N+1),'-','LineWidth',2)
            else
                plot(0:m*T-1,cell_state_dated(1+(n)*floor(T/2)-N+1:m*T+(n)*floor(T/2)-N+1),'-','LineWidth',2)
            end
        end
    end
    hold on

    % plotting cell benefit   
    nexttile(t1,2)
    if(c0_indx == 2)
    if(N <= floor(T/2))
        plot(0:m*T-1,benefit(1+floor(T/2)-N+1:m*T+floor(T/2)-N+1),'--','LineWidth',2)
    else
        n = 2;
        while(N > n*floor(T/2))
            n = n+1;
        end
        if(mod(n,2) == 0)
            plot(0:m*T-1,benefit(1+(n+1)*floor(T/2)-N+1:m*T+(n+1)*floor(T/2)-N+1),'--','LineWidth',2)
        else
            plot(0:m*T-1,benefit(1+(n)*floor(T/2)-N+1:m*T+(n)*floor(T/2)-N+1),'--','LineWidth',2)
        end
    end
    else
        if(N <= floor(T/2))
            plot(0:m*T-1,benefit(1+floor(T/2)-N+1:m*T+floor(T/2)-N+1),'-','LineWidth',2)
        else
            n = 2;
            while(N > n*floor(T/2))
                n = n+1;
            end
            if(mod(n,2) == 0)
                plot(0:m*T-1,benefit(1+(n+1)*floor(T/2)-N+1:m*T+(n+1)*floor(T/2)-N+1),'-','LineWidth',2)
            else
                plot(0:m*T-1,benefit(1+(n)*floor(T/2)-N+1:m*T+(n)*floor(T/2)-N+1),'-','LineWidth',2)
            end
        end
    end
    hold on


    nexttile(t1,3)
    q01 = (reward(2) - cost(1) + c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
    q10 = (reward(2) - cost(1) - c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
    [prob_state_1,prob_state_1_given_1] = undated_memory_cell_state_and_cond_prob_cell_state_periodic_env(N,T,q01,q10);

    plot(0:m*T-1,repmat([prob_state_1(floor(T/2)) prob_state_1(floor(T/2)+1:end) prob_state_1(1:floor(T/2)-1)],1,m),'-','LineWidth',2);
    hold on

    nexttile(t1,4)
    [undated_benefit, undated_benefit_with_time] = benefit_centered_moments_cal_periodic_env(T,reward,cost, prob_state_1, prob_state_1_given_1,c0(c0_indx)); 

    % cumulative benefit
    plot(0:m*T-1,repmat([undated_benefit_with_time(floor(T/2)+1:end); undated_benefit_with_time(1:floor(T/2))],1,m),'-','LineWidth',2);
    hold on

end

for j = 1:4
        nexttile(t1,j);
        % drawing vertical lines to compare env and cell state
        % switching (Disable this plot when considering stochastic environments)
        for i = 1:2*m-1
            plot(i*floor(T/2)*ones(10,1),linspace(0,2,10),'--k','LineWidth',1.5);
        end

        grid on
        ax = gca;
        ax.FontSize = 16;
        xlabel('time')
        if(mod(j,2)~=0)
            ylabel('P(S_n = 1)')
            ax.YTick = 0:0.25:1;
            ax.YTickLabel = 0:0.25:1;
            ylim([0 1])
        else
            ylabel('Avg Benefit')
        end
        ax.XTick = 0:2:T+1;
        xlim([0 m*T])
end
% nexttile(2);
% % drawing vertical lines to compare env and cell state
% % switching (Disable this plot when considering stochastic environments)
% for i = 1:2*m-1
%     plot(i*floor(T/2)*ones(10,1),linspace(0,1,10),'--k','LineWidth',1.5);
% end
% 
% grid on
% ax = gca;
% ax.FontSize = 16;
% xlabel('time')
% ylabel('P(S_n = 1)')
% ax.YTick = 0:0.25:1;
% ax.YTickLabel = 0:0.25:1;
% ax.XTick = 0:2:T+1;
% xlim([0 m*T])
% cd(results)

lg = legend(num2str(round(c0,2)'),'Location','northeastoutside');
title(lg,'Transition Cost');

% cd(results)
% saveas(fig1,['Dyn env and pi estimate periodic env q ' num2str(q) ' N ' num2str(N) ' T ' num2str(T) '.png'])
% saveas(fig2,['Dyn cell state cost variation periodic env q ' num2str(q) ' N ' num2str(N) ' T ' num2str(T) '.png'])

%% Cell benefit with cost variation in periodic environment

q = 0.5;
alpha = 1;
beta = alpha*q./(1-q);
cost = [1 1];
reward = [alpha+cost(2) beta+cost(1)];
N = 6;
T = 2:16;

x_high = []; % for dated memory
y_high = []; % for dated memory

x_high_1 = []; % for undated memory
y_high_1 = []; % for undated memory

w = 1;
h = 1;

c0 = linspace(0,(min(alpha,beta)-min(alpha,beta)/100),50);

benefit = zeros(length(T),length(c0),2);
% color_code= ["#eaeded","#bfc9ca","#95a5a6","#717d7e","#4d5656"];

N_indx = 1;

for T_indx = 1:length(T) 
parfor c0_indx = 1:length(c0)
    disp(['running analysis for c0 ' num2str(c0(c0_indx)) ' T ' num2str(T(T_indx))]);

    q01 = (reward(2) - cost(1) + c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
    q10 = (reward(2) - cost(1) - c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

    % % simulation of the undated memory model
    % [cell_state_undated,benefit,avg_benefit,avg_benefit_with_time,pi_est_undated] = handles{2}(reward,cost,total_time_steps, env, N, p, c0(c0_indx));

    % % simulation of the dated memory model
    % [cell_state_dated,~,avg_benefit,avg_benefit_with_time,pi_est_dated] = handles{1}(reward,cost,total_time_steps, env, N, p, c0(c0_indx));

        % determining benefit for the dated memory model
        if(N(N_indx) < T(T_indx))
            if(N(N_indx) <= floor(T(T_indx)/2))
                kmax = N(N_indx); kmin = 0;
                % for 0 -> 1 transition
                steps_for_01_trans = floor(N(N_indx)*q01)+1;
                % for 1 -> 0 transition
                steps_for_10_trans = kmax - floor(N(N_indx)*q10);
            else
                kmax = ceil(T(T_indx)/2); kmin = N(N_indx) - floor(T(T_indx)/2);
                % for 0 -> 1 transition
                steps_for_01_trans = N(N_indx) - ceil(T(T_indx)/2)+floor(N(N_indx)*q01)+1-kmin;
                % for 1 -> 0 transition
                steps_for_10_trans = N(N_indx) - floor(T(T_indx)/2) - floor(N(N_indx)*q10) + kmax;
            end
        else
            l = rem(N(N_indx),T(T_indx));
            n = (N(N_indx)-l)/T(T_indx);
            if(l <= floor(T(T_indx)/2))
                kmax = n*ceil(T(T_indx)/2)+l; kmin = n*ceil(T(T_indx)/2);
                % for 0 -> 1 transition
                steps_for_01_trans = floor(N(N_indx)*q01)+1-kmin;
                % for 1 -> 0 transition
                steps_for_10_trans = kmax - floor(N(N_indx)*q10);
            else
                kmax = (n+1)*ceil(T(T_indx)/2); kmin = n*ceil(T(T_indx)/2)+l-floor(T(T_indx)/2);
                % for 0 -> 1 transition
                steps_for_01_trans = l - ceil(T(T_indx)/2)+floor(N(N_indx)*q01)+1-kmin;
                % for 1 -> 0 transition
                steps_for_10_trans = l - floor(T(T_indx)/2) - floor(N(N_indx)*q10) + kmax;
            end
        end

        if(N(N_indx)*q01 >= kmax && N(N_indx)*q10 >= kmin)
            dated_benefit = cost(2)*floor(T(T_indx)/2)+ reward(2)*ceil(T(T_indx)/2);
        elseif(N(N_indx)*q10 < kmin && N(N_indx)*q01 < kmax)
            dated_benefit= cost(1)*floor(T(T_indx)/2)+ reward(1)*ceil(T(T_indx)/2);
        elseif(N(N_indx)*q10 < kmin && N(N_indx)*q01 >= kmax) % the asymptotic average benefit is initial cell state dependent as there won't be cell state switching 
            dated_benefit = (cost(2)*floor(T(T_indx)/2)+ reward(2)*ceil(T(T_indx)/2)); % + 0.5*(cost(1)*floor(T/2)+ reward(1)*ceil(T/2)); assuming the initial cell state is 0
        else
            dated_benefit = cost(1)*steps_for_01_trans + reward(1)*(ceil(T(T_indx)/2) - steps_for_01_trans) + cost(2)*steps_for_10_trans + reward(2)*(floor(T(T_indx)/2)-steps_for_10_trans);
        end
        dated_benefit = dated_benefit/T(T_indx);

        if(N(N_indx)*q01 >= kmax ||  N(N_indx)*q10 < kmin)
            x_high = [x_high c0_indx];
            y_high = [y_high T_indx]; % Row index - 0.5
        end

    
    [prob_state_1,prob_state_1_given_1] = undated_memory_cell_state_and_cond_prob_cell_state_periodic_env(N(N_indx),T(T_indx),q01,q10);
    undated_benefit = benefit_centered_moments_cal_periodic_env(T(T_indx),reward,cost, prob_state_1, prob_state_1_given_1,c0(c0_indx)); 
    
    benefit(T_indx,c0_indx,:) = [dated_benefit undated_benefit];
    

    centered_moments_adapt_time = undated_memory_adaptation_time_periodic_env(N,T(T_indx),q01,q10);
    if(centered_moments_adapt_time(1,1) > floor(T(T_indx)/2) || centered_moments_adapt_time(1,2) > ceil(T(T_indx)/2))
        x_high_1 = [x_high_1 c0_indx];
        y_high_1 = [y_high_1 T_indx]; % Row index - 0.5
    end

end
end

cmap = color_map_func();

figure('Position',[589 311 534 565]);
t1 = tiledlayout(2,1);
for j = 1:2
    % figure('Position',[680 629 610 249])
    nexttile
    imagesc(benefit(:,:,j),[1 2]);
    colormap(cmap)
    c = colorbar;
    % c.Limits = [0 1];
    % axis equal tight; % Adjust axis for proper grid alignment
    % Highlight specific grids (e.g., row 3, col 4)
    hold on;
    if(j == 1)
        for i = 1:length(y_high)
            rectangle('Position', [x_high(i)-0.5, y_high(i)-0.5, w, h], 'EdgeColor', 'k', 'LineWidth', 1);
        end
    else
        for i = 1:length(y_high_1)
            rectangle('Position', [x_high_1(i)-0.5, y_high_1(i)-0.5, w, h], 'EdgeColor', 'k', 'LineWidth', 1);
        end
    end
    
    % heatmap(c0,T,benefit(:,:,j),"Colormap",cmap,"ColorLimits",[1 2],'GridVisible','off');
    ax= gca;
    ax.FontSize = 16;
    ax.YTick = [1:4:length(T) length(T)];
    ax.YTickLabels = [T(1:4:length(T)) T(end)];
    ax.XTick = [1:10:50 length(c0)];
    ax.XTickLabels = [round(c0(1:10:length(c0)),2) round(c0(end),2)];
    xlabel('Transition cost')
    ylabel('Time Period (T)')
end

saveas(gcf,['Benefit variation with time period and transition cost periodic env N ' num2str(N) ' q ' num2str(q) '.png']);

% saveas(gcf,['Benefit variation with time period and transition cost periodic env N ' num2str(N) ' q ' num2str(q) ' added label undated adapt time.png']);


%% Comparison to experimental data: checking for convexity in growth; and growth reversal

close all
T = 14;
N = 8;
total_time_steps = 1000;
c0 = 0; % adaptation cost

% analysis for the undated memory model
P0 = transition_matrix_P_cal(N,0);
P1 = transition_matrix_P_cal(N,1);

% stationary distribution calculation
prob_analy = zeros(T,N+1);
prob_analy(end,:) = null((eye(N+1) - P0^(floor(T/2))*P1^(floor(T/2)+1))');
prob_analy(end,:)  = prob_analy(end,:)/sum(prob_analy(end,:));
for indx = 1:T-1
    if(indx <= floor(T/2))
        prob_analy(indx,:) = prob_analy(end,:)*P0^(indx);
    else
        prob_analy(indx,:) = prob_analy(end,:)*P0^(floor(T/2))*P1^(indx-floor(T/2));
    end
end

m = 1;
q = 0.3;%[0.3 0.5 0.7];

benefit_mat = zeros(19,2);
benefit_mat_state_1 = zeros(19,2);
benefit_mat_state_0 = zeros(19,2);
for q_indx = 1:length(q)
    T_indx = 1;
for T = 2:60
alpha = 1;
beta = alpha*q(q_indx)./(1-q(q_indx));
cost = [1 1];
reward = [alpha+cost(2) beta+cost(1)];

% color_code= ["#eaeded","#bfc9ca","#95a5a6","#717d7e","#4d5656"];

    % % simulation of the undated memory model
    % [cell_state_undated,benefit,avg_benefit,avg_benefit_with_time,pi_est_undated] = handles{2}(reward,cost,total_time_steps, env, N, p, c0(c0_indx));

    env = discrete_period_env(T,total_time_steps+N); % the state of the environment along each axis and timesteps

    % simulation of the dated memory model
    [cell_state_dated,benefit,avg_benefit,avg_benefit_with_time,pi_est_dated] = handles{1}(reward,cost,total_time_steps, env, N, p, c0);
    
    benefit_mat(T_indx,1) = avg_benefit;
    indices = find(env == 0);
    indices = indices - N;
    indices = indices(indices>0);
    benefit_mat_state_0(T_indx,1) = sum(benefit(indices))/length(indices);


    indices = find(env == 1);
    indices = indices - N;
    indices = indices(indices>0);
    benefit_mat_state_1(T_indx,1) = sum(benefit(indices))/length(indices);

    q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
    q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
    [prob_state_1,prob_state_1_given_1] = undated_memory_cell_state_and_cond_prob_cell_state_periodic_env(N,T,q01,q10);

    [undated_benefit, undated_benefit_with_time] = benefit_centered_moments_cal_periodic_env(T,reward,cost, prob_state_1, prob_state_1_given_1,c0); 
    benefit_mat(T_indx,2) = undated_benefit;
    benefit_mat_state_0(T_indx,2) = sum(undated_benefit_with_time(1:floor(T/2)))/floor(T/2);
    benefit_mat_state_1(T_indx,2) = sum(undated_benefit_with_time(floor(T/2)+1:T))/ceil(T/2);

    T_indx = T_indx+1;
end
end


% figure('Position',[680 722 540 256])
figure
plot(3:2:60,benefit_mat(2:2:59,:),'LineWidth',2);
hold on
plot(0:m*T,repmat(sum(reward)/2 * ones(1,T+1),1,m),'--','Color',"#ADA1A1",'LineWidth',2);
plot(0:m*T,repmat(reward(2)*ones(1,T+1),1,m),'--k','LineWidth',2);
plot(0:m*T,repmat(reward(1)*ones(1,T+1),1,m),'--k','LineWidth',2);
xlabel('Time period (T)');
ylabel('Average Benefit');
ax =gca;
ax.FontSize = 14;
grid on;
% saveas(gcf,['Periodic environment comparison with static environment c0 ' num2str(c0) ' N ' num2str(N) ' q ' num2str(q) '.png'])


% figure('Position',[680 722 540 256])
figure
lin_type = {'-','-.'};
for i = 1:2
plot(3:2:60,benefit_mat_state_0(2:2:59,i),'Color',"#0072BD",'LineStyle',lin_type{i},'LineWidth',2);
hold on
plot(3:2:60,benefit_mat_state_1(2:2:59,i),'Color',"#D95319",'LineStyle',lin_type{i},'LineWidth',2);
end
% plot(0:m*T-1,repmat(sum(reward)/2 * ones(1,T),1,m),'--','Color',"#ADA1A1",'LineWidth',2);
plot(0:m*T,repmat(reward(2)*ones(1,T+1),1,m),'--','Color',"#0072BD",'LineWidth',2);
plot(0:m*T,repmat(reward(1)*ones(1,T+1),1,m),'--','Color',"#D95319",'LineWidth',2);
xlabel('Time period (T)');
ylabel('Average Benefit');
ax =gca;
ax.FontSize = 14;
grid on;
% saveas(gcf,['Periodic environment component wise benefit c0 ' num2str(c0) ' N ' num2str(N) ' q ' num2str(q) '.png'])

lg = legend(num2str(q'),'Location','northeastoutside');
title(lg,'p_I');

% cd(results)
% saveas(fig1,['Dyn env and pi estimate periodic env q ' num2str(q) ' N ' num2str(N) ' T ' num2str(T) '.png'])
% saveas(fig2,['Dyn cell state indiff prob variation periodic env c0 ' num2str(c0) ' N ' num2str(N) ' T ' num2str(T) '.png'])
