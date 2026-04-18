
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                This file contains all the codes to             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate results of cellular adaptation periodic environments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Each section headline represents the analysis carried out %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Using analytical calculations to compare dynamics of cell state for specific cases
close all
N  = 11;
T = [10 20 40];
total_time_steps = 100000;
time_steps = T(end)*2;

Eh_bias = 0.5;

cost = [0.3 0.1];
reward = [0.8 1];
q = (reward(2)-cost(1))/(sum(reward) - sum(cost));
alpha = reward(1)-cost(2);
beta = alpha*q./(1-q);

c0 = 0;

q01 = (reward(2) - cost(1) + c0)/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
q10 = (reward(2) - cost(1) - c0)/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

adapt_time =zeros(length(N),length(T),2);

fig1 = figure('Position',[680 338 479 540]);
tiledlayout(3,1,"TileSpacing","compact","Padding","loose")

fig2 = figure('Position',[680 338 479 540]);
hold on
tiledlayout(3,1,"TileSpacing","compact","Padding","loose")

fig3 = figure('Position',[680 338 479 540]);
tiledlayout(3,1,"TileSpacing","compact","Padding","loose")

for N_indx = 1:length(N)

    for T_indx = 1:length(T)
        quo = floor(N(N_indx)/T(T_indx));
        n_step = (quo+1)*T(T_indx);
        
        % simulation of the undated memory model
        env = discrete_period_biased_env(T(T_indx),total_time_steps+N(N_indx),Eh_bias); % the state of the environment along each axis and timesteps
        [cell_state_undated,benefit_with_time_num,avg_benefit,avg_benefit_with_time,pi_est_undated] = Phenotypic_adaptation_models(reward,cost,total_time_steps, env, N(N_indx), 0, c0);
        

        figure(fig1);
        nexttile
        m = 1;
        num_reps = time_steps/T(T_indx);

        bar(0:num_reps*T(T_indx)-1,repmat(env(n_step+1:n_step+m*T(T_indx)),num_reps,1),0.9,'LineWidth',1.5,'FaceColor',"#77AC30",'FaceAlpha',0.7,'EdgeAlpha',0.7);
        hold on

        [prob_state_1,prob_state_1_given_1, mean_pi] = cell_state_and_cond_prob_cell_state_bias_periodic_env(N(N_indx),T(T_indx),q01,q10,Eh_bias);
        [benefit(N_indx,T_indx), benefit_with_time] = benefit_centered_moments_cal_periodic_env(T(T_indx),Eh_bias,reward,cost, prob_state_1, prob_state_1_given_1,c0);

        mean_pi = [mean_pi(floor(T(T_indx)*(1-Eh_bias))+1:end) mean_pi(1:floor(T(T_indx)*(1-Eh_bias))) ]';
        plot([0:num_reps*T(T_indx)-1],repmat(mean_pi,num_reps,1),'LineWidth',3,'Color',"#0072BD");
        hold on
        plot([0:num_reps*T(T_indx)-1],repmat(q,time_steps,1),'--','LineWidth',3,'Color','r');


        figure(fig2)
        nexttile
        anal_cell_state_prob = [prob_state_1(floor(T(T_indx)*(1-Eh_bias))) prob_state_1(floor(T(T_indx)*(1-Eh_bias))+1:end) prob_state_1(1:floor(T(T_indx)*(1-Eh_bias))-1)];

        plot([0:num_reps*T(T_indx)-1],repmat(anal_cell_state_prob,1,num_reps),'LineWidth',3,'Color',"#0072BD");
        hold on

        figure(fig3);
        benefit_with_time_anal = [benefit_with_time(floor(T(T_indx)*(1-Eh_bias))+1:end); benefit_with_time(1:floor(T(T_indx)*(1-Eh_bias)))];

        nexttile
        plot([0:num_reps*T(T_indx)-1],repmat(benefit_with_time_anal,num_reps,1),'LineWidth',3,'Color',"#0072BD");
        hold on
        plot([0:num_reps*T(T_indx)-1],repmat(benefit(N_indx,T_indx),time_steps,1),'--','LineWidth',3,'Color','r');

        m = 1000;
        n_step = (quo+1)*T(T_indx) - N(N_indx);

        cell_state_undated =  cell_state_undated(n_step+1:n_step+m*T(T_indx));
        mean_cell_state_undated = mean(reshape(cell_state_undated',T(T_indx),[]),2);

        benefit_with_time_num = benefit_with_time_num(n_step+1:n_step+m*T(T_indx));
        mean_benefit_with_time_num = mean(reshape(benefit_with_time_num',T(T_indx),[]),2);

        figure(fig2);
        plot(0:num_reps*T(T_indx)-1,repmat(mean_cell_state_undated,num_reps,1),'o','MarkerEdgeColor',	"#D95319",'MarkerSize',8,'LineWidth',1)

        figure(fig3)
        plot(0:num_reps*T(T_indx)-1,repmat(mean_benefit_with_time_num,num_reps,1),'o','MarkerEdgeColor',	"#D95319",'MarkerSize',8,'LineWidth',1)

        figure(fig1);
        hold on
        grid on
        ax = gca;
        ax.FontSize = 16;
        if(T_indx==length(T))
            xlabel('time')
        end
        ylabel('P(\pi_n)')
        ax.XTick = 0:floor(T(1)):time_steps-1;

        figure(fig2);
        hold on
        grid on
        ax = gca;
        ax.FontSize = 16;
        if(T_indx==length(T))
            xlabel('time')
        end
        ylabel('P(S_n)')

        ax.XTick = 0:floor(T(1)):time_steps-1;
        ylim([0 1])
        
        figure(fig3);
        hold on

        grid on
        ax = gca;
        ax.FontSize = 16;
        if(T_indx==length(T))
            xlabel('time')
        end
        ylabel('E(R_n)')
        ax.XTick = 0:floor(T(1)):time_steps-1;
        ylim([min(cost) max(reward)])
    end
end
% saveas(fig1,['periodic_cell_pi_est_N_' num2str(N) '_T_10_20_40_Eh_bias_' num2str(Eh_bias) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png']);
% saveas(fig2,['periodic_cell_state_N_' num2str(N) '_T_10_20_40_Eh_bias_' num2str(Eh_bias) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png']);
% saveas(fig3,['periodic_cell_benefit_N_' num2str(N) '_T_10_20_40_Eh_bias_' num2str(Eh_bias) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png']);

%% mean benefit in periodic environments
N  = 11:11:33;
T = 20:120;
total_time_steps = 100000;
Eh_bias = 0.5;

cost = [0.3 0.1];
reward = [0.8 1];
q = (reward(2)-cost(1))/(sum(reward) - sum(cost));
alpha = reward(1)-cost(2);
beta = alpha*q./(1-q);

c0 = linspace(0,(min(alpha,beta)-min(alpha,beta)/2),3);

benefit = zeros(length(N),length(c0),length(T),length(Eh_bias),2);
adaptation_time = zeros(2,2,length(N),length(c0),length(T),length(Eh_bias));

Eh_bias_len = length(Eh_bias);
c0_len = length(c0);
T_len = length(T);
for N_indx = 1:length(N)
    for T_indx = 1:T_len
        for bias_indx = 1:Eh_bias_len
            for c0_indx = 1:c0_len
                q01 = (reward(2) - cost(1) + c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
                q10 = (reward(2) - cost(1) - c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
                
                env = discrete_period_biased_env(T(T_indx),total_time_steps+N(N_indx),Eh_bias(bias_indx)); % the state of the environment along each axis and timesteps
                [cell_state_undated,~,avg_benefit,avg_benefit_with_time,pi_est_undated] = Phenotypic_adaptation_models(reward,cost,total_time_steps, env, N(N_indx), 0, c0(c0_indx));
                
                [prob_state_1,prob_state_1_given_1] = cell_state_and_cond_prob_cell_state_bias_periodic_env(N(N_indx),T(T_indx),q01,q10,Eh_bias(bias_indx));
                benefit(N_indx,c0_indx,T_indx,bias_indx,:) = [avg_benefit benefit_centered_moments_cal_periodic_env(T(T_indx),Eh_bias(bias_indx),reward,cost, prob_state_1, prob_state_1_given_1,c0(c0_indx))];
            end
        end
    end
end

%% benefit variation with mem size

figure('Position',[680 643 438 235])
color_code = ["#0072BD", "#D95319", "#77AC30"];
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
c0_indx = 1;
for bias_indx = 1%[3 5 7]%1:length(Eh_bias)

    nexttile
    plt_indx = 1;
    for T_indx = 1:length(T)

        for m = 1:2
            if (m==1)

                plot(N(1:2:end),reshape(benefit(1:2:end,c0_indx,T_indx,bias_indx,m),length(N(1:2:end)),1),'o','MarkerEdgeColor',color_code(plt_indx),'MarkerSize',8,'LineWidth',1);
                hold on
            else
                plot(N,reshape(benefit(:,c0_indx,T_indx,bias_indx,m),length(N),1),'LineWidth',2,'Color',color_code(plt_indx));
                hold on
            end
        end
        plt_indx = plt_indx+1;
    end
    grid on
    ax =gca;
    ax.FontSize = 14;
    ax.XTick = 2:4:N(end);
    if(bias_indx == 1)
        xlabel('memory size (m)');
    end
    ylabel('Avg Benefit');
    xlim([N(1) N(end)])
end
 

% saveas(gcf,['periodic_cell_benefit_variation_with_mem_size_T_10_20_40_Eh_bias_' num2str(Eh_bias) ' match reward ' num2str(1./reward) ' mismatch reward ' num2str(1./cost) '_c0_' num2str(c0) '.png']);

%% benefit variation with adaptation cost

figure('Position',[680 643 438 235])
color_code = ["#0072BD", "#D95319", "#77AC30"];
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
N_indx = 1;
for bias_indx = 1:length(Eh_bias)
    nexttile
    plt_indx = 1;
    for T_indx = 1:length(T)

        for m = 1:2
            if (m==1)
                    plot(c0(1:2:end),reshape(benefit(N_indx,1:2:end,T_indx,bias_indx,m),length(c0(1:2:end)),1),'o','MarkerEdgeColor',color_code(plt_indx),'MarkerSize',8,'LineWidth',1);
                    hold on
            else
                plot(c0,reshape(benefit(N_indx,:,T_indx,bias_indx,m),length(c0),1),'LineWidth',2,'Color',color_code(plt_indx));
                hold on
            end
        end
        plt_indx = plt_indx+1;
    end
    grid on
    ax =gca;
    ax.FontSize = 14;
    % ax.XTick = 1:1:N(end);
    if(bias_indx == 1)
        xlabel('adaptation cost');
    end
    ylabel('Avg fitness');
    xlim([c0(1) c0(end)])
end
% saveas(gcf,['periodic_cell_benefit_variation_with_cost_T_10_20_40_Eh_bias_' num2str(Eh_bias) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '_N_' num2str(N) '.png']);

%% benefit variation with time period

color_code = ["#0072BD", "#D95319", "#77AC30"];

for c0_indx = 1:length(c0)
    figure('Position',[680 643 438 235])
    tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
for bias_indx = 1%[3 5 7]%1:length(Eh_bias)

    nexttile
    plt_indx = 1;
    for N_indx = 1:length(N)

        for m = 1:2
            if (m==1)

                % plot(T(1:2:end),reshape(benefit(N_indx,c0_indx,1:2:end,bias_indx,m),length(T(1:2:end)),1),'o','MarkerEdgeColor',color_code(plt_indx),'MarkerSize',8,'LineWidth',1);
                hold on
            else
                plot(T,reshape(benefit(N_indx,c0_indx,:,bias_indx,m),length(T),1),'LineWidth',2,'Color',color_code(plt_indx));
                hold on
            end
        end
        plt_indx = plt_indx+1;
    end
    grid on
    ax =gca;
    ax.FontSize = 14;
    ax.XTick = T(1):10:T(end);
    if(bias_indx == 1)
        xlabel('time period (T)');
    end
    ylabel('Avg fitness');
    xlim([T(1) T(end)])
end
% saveas(gcf,['periodic_cell_benefit_variation_with_time_period_N_11_22_33_Eh_bias_' num2str(Eh_bias) ' match reward ' num2str(1./reward) ' mismatch reward ' num2str(1./cost) '_c0_' num2str(c0(c0_indx)) '.png']);
end

%% residence time variation with memory size
color_code = ["#0072BD", "#D95319", "#77AC30"];
c0_indx = 1;
state_label = ['l','h'];
for cell_state = 1:2
    figure('Position',[680 643 438 235])
    tiledlayout(1,1,"TileSpacing","compact","Padding","compact")

    for bias_indx = 1%[3 5 7]%1:length(Eh_bias)
        nexttile
        plt_indx = 1;
        for T_indx = 1:length(T)

            for m = 2%1:2
                if (m==1)
                    % plot(N(1:2:end),reshape(benefit(1:2:end,c0_indx,T_indx,bias_indx,m),length(N(1:2:end)),1),'o','MarkerEdgeColor',color_code(plt_indx),'MarkerSize',8,'LineWidth',1);
                    hold on
                else
                    plot(N,reshape(adaptation_time(1,cell_state,:,c0_indx,T_indx,bias_indx),length(N),1),'LineWidth',2,'Color',color_code(plt_indx));
                    hold on
                end
            end
            plt_indx = plt_indx+1;
        end
        grid on
        ax =gca;
        ax.FontSize = 14;
        ax.XTick = 2:4:N(end);
        if(bias_indx == 1)
            xlabel('memory size (m)');
        end
        ylabel([{'Mean Residence'},{['State S^' state_label(cell_state)]}])
        xlim([N(1) N(end)])
    end
    saveas(gcf,['residence_time_variation_with_mem_size_periodic_env_state_S_' num2str(state_label(cell_state)) '_c0_' num2str(c0_set(1)) '_Eh_bias_' num2str(Eh_bias) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
end

%% residence time variation with adaptation cost
color_code = ["#0072BD", "#D95319", "#77AC30"];
N_indx = 1;
state_label = ['l','h'];
for cell_state = 1:2
    figure('Position',[680 643 438 235])
    tiledlayout(1,1,"TileSpacing","compact","Padding","compact")

    for bias_indx = 1%[3 5 7]%1:length(Eh_bias)
        nexttile
        plt_indx = 1;
        for T_indx = 1:length(T)

            for m = 2%1:2
                if (m==1)
                    % plot(N(1:2:end),reshape(benefit(1:2:end,c0_indx,T_indx,bias_indx,m),length(N(1:2:end)),1),'o','MarkerEdgeColor',color_code(plt_indx),'MarkerSize',8,'LineWidth',1);
                    hold on
                else
                    plot(c0,reshape(adaptation_time(1,cell_state,N_indx,:,T_indx,bias_indx),length(c0),1),'LineWidth',2,'Color',color_code(plt_indx));
                    hold on
                end
            end
            plt_indx = plt_indx+1;
        end
        grid on
        ax =gca;
        ax.FontSize = 14;
        % ax.XTick = 2:4:c0(end);
        if(bias_indx == 1)
            xlabel('adaptation cost');
        end
        ylabel([{'Mean Residence'},{['State S^' state_label(cell_state)]}])
        xlim([c0(1) c0(end)])
        yscale log
    end
    saveas(gcf,['residence_time_variation_with_adaptation_cost_periodic_env_state_S_' num2str(state_label(cell_state)) '_N_' num2str(N) '_Eh_bias_' num2str(Eh_bias) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
end

%% time dynamics of cell state with mem size variation

m = 1; % number of time period to be plotted
c0 = 0;%linspace(0,(min(alpha,beta)-min(alpha,beta)/100),5);
c0_indx = 1;
figure('Position',[680 706 752 172]);
t1 = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
Eh_bias = 0.5;
cost = [0.3 0.1];
reward = [0.8 1];
q = (reward(2)-cost(1))/(sum(reward) - sum(cost));
alpha = reward(1)-cost(2);
beta = alpha*q./(1-q);

reward = [alpha+cost(2) beta+cost(1)];

T = 20;
N = [5 10 20 40];

for N_indx = 1:length(N)
    nexttile(t1,1)
    q01 = (reward(2) - cost(1) + c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
    q10 = (reward(2) - cost(1) - c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
    [prob_state_1,prob_state_1_given_1,mean_pi] = cell_state_and_cond_prob_cell_state_bias_periodic_env(N(N_indx),T,q01,q10,Eh_bias);

    plot(0:m*T-1,repmat([prob_state_1(floor(T*(1-Eh_bias))) prob_state_1(floor(T*(1-Eh_bias))+1:end) prob_state_1(1:floor(T*(1-Eh_bias))-1)],1,m),'-','LineWidth',2);
    hold on

    nexttile(t1,2)
    [undated_benefit, benefit_with_time] = benefit_centered_moments_cal_periodic_env(T,Eh_bias,reward,cost, prob_state_1, prob_state_1_given_1,c0(c0_indx)); 

    % cumulative benefit
    plot(0:m*T-1,repmat([benefit_with_time(floor(T*(1-Eh_bias))+1:end); benefit_with_time(1:floor(T*(1-Eh_bias)))],1,m),'-','LineWidth',2);
    hold on

end

for j = 1:2
        nexttile(t1,j);
        for i = 1:2*m-1
            plot(i*floor(T*(Eh_bias))*ones(10,1),linspace(-4,3,10),'--k','LineWidth',1.5);
        end

        grid on
        ax = gca;
        ax.FontSize = 10;
        xlabel('time')
        if(j == 1)
            ylabel('P(S_n = S^h)')
            ax.YTick = 0:0.25:1;
            ax.YTickLabel = 0:0.25:1;
            ylim([0 1])
        else
            ylabel('Avg Benefit')
            ylim([0 1])
            ax.YTick = (0:0.25:1);
        end
        % ax.FontSize = 12;
        ax.XTick = 0:2:T+1;
        xlim([0 m*T])
end

lg = legend(num2str(N'),'Location','northeastoutside');
title(lg,'memory size (m)');
% saveas(gcf,['avg_benefit_cell_state_with_time_mem_size_variation_T_' num2str(T) '_Eh_bias_' num2str(Eh_bias) '_c0_' num2str(c0) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
%% time dynamics of cell state with cost variation

m = 1; % number of time period to be plotted
figure('Position',[680 706 752 172]);
t1 = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

Eh_bias = 0.5;

cost = [0.3 0.1];
reward = [0.8 1];
q = (reward(2)-cost(1))/(sum(reward) - sum(cost));
alpha = reward(1)-cost(2);
beta = alpha*q./(1-q);
c0 = [0 0.1 0.3 0.4 0.5];%linspace(0,(min(alpha,beta)-min(alpha,beta)/100),5);

T = 20;
N = 11;
for c0_indx = 1:length(c0)

    nexttile(t1,1)
    q01 = (reward(2) - cost(1) + c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
    q10 = (reward(2) - cost(1) - c0(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition
    [prob_state_1,prob_state_1_given_1,mean_pi] = cell_state_and_cond_prob_cell_state_bias_periodic_env(N,T,q01,q10,Eh_bias);

    plot(0:m*T-1,repmat([prob_state_1(floor(T*(1-Eh_bias))) prob_state_1(floor(T*(1-Eh_bias))+1:end) prob_state_1(1:floor(T*(1-Eh_bias))-1)],1,m),'-','LineWidth',2);
    hold on

    nexttile(t1,2)
    [undated_benefit, benefit_with_time] = benefit_centered_moments_cal_periodic_env(T,Eh_bias,reward,cost, prob_state_1, prob_state_1_given_1,c0(c0_indx)); 

    % cumulative benefit
    plot(0:m*T-1,repmat([benefit_with_time(floor(T*(1-Eh_bias))+1:end); benefit_with_time(1:floor(T*(1-Eh_bias)))],1,m),'-','LineWidth',2);
    hold on

end

for j = 1:2
        nexttile(t1,j);
        for i = 1:2*m-1
            plot(i*floor(T*(Eh_bias))*ones(10,1),linspace(0,2,10),'--k','LineWidth',1.5);
        end

        grid on
        ax = gca;
        ax.FontSize = 10;
        xlabel('time')
        if(j == 1)
            ylabel('P(S_n = S^h)')
            ax.YTick = (0:0.25:1);
            ax.YTickLabel = (0:0.25:1);
            ylim([0 1])
        else
            ylabel('Avg Benefit')
            ylim([0 1])
            ax.YTick = (0:0.25:1);
        end
        % ax.FontSize = 12;
        ax.XTick = 0:2:T+1;
        xlim([0 m*T])
end

lg = legend(num2str(round(c0,2)'),'Location','northeastoutside');
title(lg,'Adaptation Cost');

% saveas(gcf,['avg_benefit_cell_state_with_time_cost_variation_T_' num2str(T) '_Eh_bias_' num2str(Eh_bias) '_N_' num2str(N) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])