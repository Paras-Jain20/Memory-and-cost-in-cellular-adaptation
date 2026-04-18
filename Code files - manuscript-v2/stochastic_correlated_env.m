% close all

total_time_steps = 1000000;

state_label = {'l', 'h'};
model_label = {'Undated', 'Dated'};

N = 2:12; % memory capacity
p = 0.5; % stochastic environmental state
num_rep= 3; % number of replicate for numberical simulations

cost = [0.3 0.1];
reward = [0.8 1];
% cost = [0 -1];
% reward = [0.2 1];
q = (reward(2)-cost(1))/(sum(reward) - sum(cost));
alpha = reward(1)-cost(2);
beta = alpha*q./(1-q);

% q = 0.55657; % indifference environment 
% 
% alpha = 1; % difference between (r_h^+)-(r_l^-) 
% beta = alpha*q./(1-q); % (r_l^+)-(r_h^-)
% cost = [0.2 0.8]; % (r_h^-), (r_l^-)
% reward = [alpha+cost(2) beta+cost(1)]; % (r_h^+), (r_l^+)

c0_set = 0;%linspace(0,min(alpha,beta)-min(alpha,beta)/100,50); % adaptation cost
% c0_set(end)
% one_step_auto_corr = linspace(-0.9,0.9,37);
one_step_auto_corr = linspace(-0.6,0.6,3);

%% Analytical calculations

transition_syn = zeros(2,2,length(c0_set),length(one_step_auto_corr),length(N));
resi_time_centered_moments = zeros(4,2,length(c0_set),length(one_step_auto_corr),length(N)); % last two indices corresponds to the number of p and model values
benefit_centered_moments = zeros(2,length(c0_set),length(one_step_auto_corr),length(N)); % last two indices corresponds to the number of p and model values
P_Sn_En = zeros(4,length(c0_set),length(one_step_auto_corr),length(N)); % last two indices corresponds to the number of p and model values
env_cell_state_corr = zeros(length(c0_set),length(one_step_auto_corr),length(N));
joint_prob_residence_time_state_moments = zeros(length(c0_set),6); % first four columns represent Expectation and Variance of x/(x+y) and y/(x+y); next two columns represent E(x+y) and Var(x+y)

prob_residence_time_state_0 = cell(length(c0_set));
prob_residence_time_state_1 = cell(length(c0_set));

joint_prob_flag = zeros(2,length(c0_set));

for N_indx = 1:length(N)
    for corr_env_indx = 1:length(one_step_auto_corr)
        p2 = min(1,(1-one_step_auto_corr(corr_env_indx))*(1-p));
        p1 = min(1,(1-one_step_auto_corr(corr_env_indx))*(p));
        P = transition_matrix_P_cal_corr_env(N(N_indx),p1,p2);

        for c0_indx = 1:length(c0_set)
            q01 = (reward(2) - cost(1) + c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
            q10 = (reward(2) - cost(1) - c0_set(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

            disp(['running analytical calculations for N ' num2str(N(N_indx)) ' c0 ' num2str(c0_set(c0_indx)) ' one step autocorr ' num2str(one_step_auto_corr(corr_env_indx)) ' independent p ' num2str(p1/(p1+p2))]);

            [resi_time_centered_moments(:,:,c0_indx,corr_env_indx,N_indx), cond_prob_Sn_Sn_1, x_1_ss,x_0_ss, x_Sn_Sn_1, P_s] = marginal_resi_time_and_cond_prob_cell_state_stoc_env(N(N_indx),q01,q10,P);
            transition_syn(:,:,c0_indx,corr_env_indx,N_indx) = [P_s(1,2) p1; P_s(2,1) p2];
            [benefit_centered_moments(:,c0_indx,corr_env_indx,N_indx),env_cell_state_corr(c0_indx,corr_env_indx,N_indx), P_Sn_En(:,c0_indx,corr_env_indx,N_indx)] = benefit_centered_moments_cal(p1,p2,reward, cost,cond_prob_Sn_Sn_1, x_1_ss,x_0_ss, x_Sn_Sn_1,c0_set(c0_indx));        
        end
    end
end
%% Numerical Calculations
c0_set_num = 0;%linspace(0,min(alpha,beta)-min(alpha,beta)/100,20);

num_benefit_centered_moments = zeros(4,length(c0_set_num),length(one_step_auto_corr),length(N),num_rep); % first two indices corresponds to the number of p and model values
num_resi_time_centered_moments = zeros(4,2,length(c0_set_num),length(one_step_auto_corr),length(N),num_rep); % last two indices corresponds to the number of p and model values

pi_est = zeros(total_time_steps,length(one_step_auto_corr),num_rep);
env_mat = zeros(total_time_steps,length(one_step_auto_corr),num_rep);

N_length = length(N);
c0_length = length(c0_set_num);
one_step_auto_corr_len = length(one_step_auto_corr);

for rep_indx = 1:num_rep
    for N_indx = 1:N_length
        temp_num_benefit_centered_moments = zeros(4,length(c0_set_num),length(one_step_auto_corr));
        temp_num_resi_time_centered_moments = zeros(4,2,length(c0_set_num),length(one_step_auto_corr));

        for c0_indx = 1:c0_length
            for corr_env_indx = 1:one_step_auto_corr_len

                q01 = (reward(2) - cost(1) + c0_set_num(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 0 -> 1 cell state transition
                q10 = (reward(2) - cost(1) - c0_set_num(c0_indx))/(sum(reward) - sum(cost)); % indifference probability for 1 -> 0 cell state transition

                temp_cell_state = zeros(total_time_steps,num_rep);
                temp_pi_est = zeros(total_time_steps,num_rep);
                disp(['running numerical analyses for N ' num2str(N(N_indx)) ' c0 ' num2str(c0_set_num(c0_indx)) ' one step autocorr ' num2str(one_step_auto_corr(corr_env_indx)) ' rep ' num2str(rep_indx)])

                p2 = min(1,(1-one_step_auto_corr(corr_env_indx))*(1-p));
                p1 = min(1,(1-one_step_auto_corr(corr_env_indx))*(p));

                env = correlated_env(p1,p2,total_time_steps+N(N_indx)); % the state of the environment along each axis and timesteps

                [temp_cell_state(:,rep_indx),benefit,avg_benefit,avg_benefit_with_time,temp_pi_est(:,rep_indx)] = Phenotypic_adaptation_models(reward,cost,total_time_steps, env, N(N_indx), p, c0_set_num(c0_indx));
                temp_auto_corr = autocorr(benefit);
                temp_num_benefit_centered_moments(:,c0_indx,corr_env_indx) = [avg_benefit; var(benefit); skewness(benefit); temp_auto_corr(2)]; % saving lag 1 autocorrelation

                [residence_time_state_0,residence_time_state_1] = residence_time_cal_cell_state_based(temp_cell_state(:,rep_indx));
                if(length(residence_time_state_0) >= 3 || length(residence_time_state_1) >= 3)
                    temp_num_resi_time_centered_moments(:,:,c0_indx,corr_env_indx) = [mean(residence_time_state_0) mean(residence_time_state_1); var(residence_time_state_0) var(residence_time_state_1); skewness(residence_time_state_0) skewness(residence_time_state_1); kurtosis(residence_time_state_0) kurtosis(residence_time_state_1)];
                else
                    temp_num_resi_time_centered_moments(:,:,c0_indx,corr_env_indx) = zeros(4,2);
                end
                pi_est(:,corr_env_indx,rep_indx) = temp_pi_est(:,rep_indx);
                env_mat(:,corr_env_indx,rep_indx) = env(N(N_indx)+1:end);
            end
        end
        temp_size = size(num_benefit_centered_moments);

        num_benefit_centered_moments(:,:,:,N_indx,rep_indx) = temp_num_benefit_centered_moments;
        num_resi_time_centered_moments(:,:,:,:,N_indx,rep_indx) = temp_num_resi_time_centered_moments;
    end

end

%% Comparing the distribution of environmental inference between two memories for increasing env correlation
figure('Position',[646 299 759 545])
tiledlayout(3,4,"TileSpacing","compact","Padding","loose");
N_indx = 1;
k = 0:N(N_indx);
for corr_env_indx = 1:3

    nexttile(4*(corr_env_indx-1)+1,[1 3])
    bar(1:50,env_mat(1:50,corr_env_indx,1,1),'FaceColor',"#77AC30")
    hold on
    plot(1:50,pi_est(1:50,corr_env_indx,1,1),'LineWidth',1.5)
    ax = gca;
    ax.FontSize = 14;
    if(corr_env_indx ~= 3)
        ax.XTickLabel = " ";
    end

    if(corr_env_indx == 3)
        xlabel('time')
    end
    ylabel('\pi_n')
    ylim([-0.05 1.05])
    ax.YTick = 0:.25:1;
    grid on
    %%%%%%%%%%%%%%%%%%%%%%
    nexttile(4+4*(corr_env_indx-1))

    p2 = min(1,(1-one_step_auto_corr(corr_env_indx))*(1-p));
    p1 = min(1,(1-one_step_auto_corr(corr_env_indx))*(p));
    P = transition_matrix_P_cal_corr_env(N(N_indx),p1,p2);
    prob_k = [1 zeros(1, size(P,1)-1)]*P^1000;
    prob_k = sum(reshape(prob_k,2,[]),1);

    bar(k/N(N_indx),prob_k,'FaceAlpha',0.99,'EdgeColor','none','Horizontal','on');

    hold on

    k_prob = zeros(N(N_indx)+1,1);
    for indx = 1:total_time_steps
        for temp_k = 0:N(N_indx)
            if(pi_est(indx,corr_env_indx,1,1)*N(N_indx) == temp_k)
                k_prob(temp_k+1) = k_prob(temp_k+1)+1;
            end
        end
    end
    k_prob = k_prob/total_time_steps;

    plot(k_prob,(0:N(N_indx))/N(N_indx),'Marker','*','MarkerSize',8,'Color','k','LineStyle','none')

    ax = gca;
    ax.FontSize = 14;
    % ax.YTick = k/N(N_indx);
    ax.YTick = 0:.25:1;
    ax.YTickLabel = " ";
    if(corr_env_indx ~= 3)
        ax.XTickLabel = " ";
    end
    if(corr_env_indx == 3)
        xlabel('P(\pi_n)')
    end
    grid on
    hold on
    ylim([-0.05 1.05])
    ax.XTick = 0:0.1:0.4;
    xlim([0 0.4])

end

%% Benefit variation with memmory size 
color_code = ["#0072BD","#D95319", "#77AC30"];
Marker_type = {'o','d'};

figure('Position',[680 550 397 232.6000])
c0_indx = 1;
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
max_y = max(benefit_centered_moments(2,c0_indx,:,:),[],'all');
min_y = min(benefit_centered_moments(2,c0_indx,:,:),[],'all');

for corr_indx = 1:3
    data = reshape(benefit_centered_moments(2,c0_indx,corr_indx,:),1,length(N));
    plot(N,data,'LineWidth',2,'Color',color_code(corr_indx));
    hold on;

    data = reshape(num_benefit_centered_moments(2,c0_indx,corr_indx,:,1),1,length(N));
    plot(N(1:end),data(1:end),'Marker',Marker_type{1},'Color',color_code(corr_indx),'MarkerSize',10,'LineStyle','none');
    
    % ylabel('Avg fitness')
    ylabel('Var fitness')

    grid on
    % ylabel('\rho_{xx}(1)')
    ax = gca;
    ax.FontSize = 14;
    ylim([min_y-0.01 max_y+0.01])
    xlim([min(N) max(N)])
    ax.YTick = round(min_y,2):(round(max_y,2)-round(min_y,2))/5:round(max_y,2);
end

xlabel('memory size (m)')
saveas(gcf,['var_benefit_variation_with_mem_size_stoc_corr_env_c0_' num2str(c0_set) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])

% saveas(gcf,['avg_benefit_variation_with_mem_size_stoc_corr_env_c0_' num2str(c0_set) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])

%% Benefit variation with cost
color_code = ["#0072BD","#D95319", "#77AC30"];

Marker_type = {'o','d','+'};

N_indx = 1;
for out_type = 1:2
figure('Position',[680 550 397 232.6000])
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
max_y = max(benefit_centered_moments(out_type,:,:,N_indx,:),[],'all');
min_y = min(benefit_centered_moments(out_type,:,:,N_indx,:),[],'all');

for corr_indx = 1:3
    data = reshape(benefit_centered_moments(out_type,:,corr_indx,N_indx),1,length(c0_set));
    plot(c0_set,data,'LineWidth',2,'Color',color_code(corr_indx));
    hold on;

    data = reshape(num_benefit_centered_moments(out_type,:,corr_indx,N_indx,1),1,length(c0_set_num));
    plot(c0_set_num,data,'Marker',Marker_type{1},'Color',color_code(corr_indx),'MarkerSize',10,'LineStyle','none');

    if(out_type == 1)
            ylabel('Avg fitness')
    else    
        ylabel('Var fitness')
    end
    grid on
    ax = gca;
    ax.FontSize = 14;
    ylim([min_y-0.01 max_y+0.01])
    xlim([min(c0_set) max(c0_set)])
    % ax.YTick = round(min_y,2):(round(max_y,2)-round(min_y,2))/5:round(max_y,2);
end

xlabel('adaptation cost')
if(out_type == 1)
saveas(gcf,['avg_benefit_variation_with_cost_stoc_corr_env_N_' num2str(N) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
else
    saveas(gcf,['var_benefit_variation_with_cost_stoc_corr_env_N_' num2str(N) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
end
end
%% Residence time variation with memmory size 
color_code = ["#0072BD","#D95319", "#77AC30"];
Marker_type = {'o','d'};
state_label = ['l','h'];

for cell_state = 1:2
figure('Position',[680 550 397 232.6000])
c0_indx = 1;
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
max_y = max(resi_time_centered_moments(1,cell_state,c0_indx,:,:),[],'all');
min_y = min(resi_time_centered_moments(1,cell_state,c0_indx,:,:),[],'all');

for corr_indx = 1:3
    data = reshape(resi_time_centered_moments(1,cell_state,c0_indx,corr_indx,:),1,length(N));
    plot(N,data,'LineWidth',2,'Color',color_code(corr_indx));
    hold on;

    data = reshape(mean(num_resi_time_centered_moments(1,cell_state,c0_indx,corr_indx,:,1),6),1,length(N));
    plot(N(1:end),data(1:end),'Marker',Marker_type{1},'Color',color_code(corr_indx),'MarkerSize',10,'LineStyle','none');
    
    % ylabel('Avg fitness')
    ylabel([{'Mean Residence'},{['State S^' state_label(cell_state)]}])

    grid on
    % ylabel('\rho_{xx}(1)')
    ax = gca;
    ax.FontSize = 14;
    ylim([min_y-1 max_y+1])
    xlim([min(N) max(N)])
    % ax.YTick = round(min_y):(round(max_y)-round(min_y,2))/5:round(max_y,2);
end

xlabel('memory size (m)')
% saveas(gcf,['residence_time_variation_with_mem_size_stoc_corr_env_state_S_' num2str(state_label(cell_state)) '_c0_' num2str(c0_set) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
end
%% Residence time variation with adaptation cost 
color_code = ["#0072BD","#D95319", "#77AC30"];
Marker_type = {'o','d'};
state_label = ['l','h'];

for cell_state = 1:2
figure('Position',[680 550 397 232.6000])
N_indx = 1;
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
max_y = max(resi_time_centered_moments(1,cell_state,:,:,N_indx),[],'all');
min_y = min(resi_time_centered_moments(1,cell_state,:,:,N_indx),[],'all');

for corr_indx = 1:3
    data = reshape(resi_time_centered_moments(1,cell_state,:,corr_indx,N_indx),1,length(c0_set));
    plot(c0_set,data,'LineWidth',2,'Color',color_code(corr_indx));
    hold on;

    data = reshape(mean(num_resi_time_centered_moments(1,cell_state,:,corr_indx,N_indx,1),6),1,length(c0_set_num));
    plot(c0_set_num(1:end),data(1:end),'Marker',Marker_type{1},'Color',color_code(corr_indx),'MarkerSize',10,'LineStyle','none');
    
    % ylabel('Avg fitness')
    ylabel([{'Mean Residence'},{['State S^' state_label(cell_state)]}])

    grid on
    % ylabel('\rho_{xx}(1)')
    ax = gca;
    ax.FontSize = 14;
    ylim([min_y-1 max_y+1])
    xlim([min(c0_set) max(c0_set)])
    yscale log
    % ax.YTick = round(min_y):(round(max_y)-round(min_y,2))/5:round(max_y,2);
end

xlabel('adaptation cost')
saveas(gcf,['residence_time_variation_with_adaptation_cost_stoc_corr_env_state_S_' num2str(state_label(cell_state)) '_N_' num2str(N) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
end

%% Generating data frame
data_table = zeros(length(c0_set)*length(N)*length(one_step_auto_corr),12); % 2 for number of models
for N_indx = 1:length(N)
    for c0_indx = 1:length(c0_set)
        for corr_env_indx = 1:length(one_step_auto_corr)
            data_table(corr_env_indx+(c0_indx-1)*length(one_step_auto_corr) + (N_indx-1)*length(one_step_auto_corr)*length(c0_set),:) = [N(N_indx) c0_indx (corr_env_indx) ...
                benefit_centered_moments(1:2,c0_indx,corr_env_indx,N_indx)' resi_time_centered_moments(1,:,c0_indx,corr_env_indx,N_indx) transition_syn(1,:,c0_indx,corr_env_indx,N_indx) ...
                transition_syn(2,:,c0_indx,corr_env_indx,N_indx) env_cell_state_corr(c0_indx,corr_env_indx,N_indx)];
        end
    end
end


data_table = splitvars(table(data_table));
% data_table.Properties.VariableNames = {'model','N','c0','auto_corr','avg_benefit','var_benefit','skew_benefit','auto_benefit','Resi_0','Resi_1'};
data_table.Properties.VariableNames = {'N','c0_indx','auto_corr_indx','avg_benefit','var_benefit','Resi_0','Resi_1','t01','p1','t10','p2','env_cell_corr'};
data_table{isnan(data_table.avg_benefit),4:12} = nan;

%% distribution of avg benefit across memory size and adaptation cost for increasing environment correlation
data_filtered = data_table(logical((sum(data_table.auto_corr_indx == [1:2:length(one_step_auto_corr)],2))),:);

for out_type = 1:2 % 1: avg benefit; 2: var benefit


temp_benefit = zeros(length(one_step_auto_corr(1:2:end)),4,2);
indx = 1;
for auto_corr_indx = 1:2:length(one_step_auto_corr)
    temp_benefit(indx,1,:) = data_filtered{logical((data_filtered.auto_corr_indx == auto_corr_indx) .* (data_filtered.N == 2) .* (data_filtered.c0_indx  == 1)), 4:5};
    temp_benefit(indx,2,:) = data_filtered{logical((data_filtered.auto_corr_indx  == auto_corr_indx) .* (data_filtered.N == N(end)) .* (data_filtered.c0_indx  == 1)), 4:5};
    temp_benefit(indx,3,:) = data_filtered{logical((data_filtered.auto_corr_indx  == auto_corr_indx) .* (data_filtered.N == 2) .* (data_filtered.c0_indx  == length(c0_set))), 4:5};
    temp_benefit(indx,4,:) = data_filtered{logical((data_filtered.auto_corr_indx  == auto_corr_indx) .* (data_filtered.N == N(end)) .* (data_filtered.c0_indx  == length(c0_set))), 4:5};

    indx = indx+1;
end

marker_type = {'<','*','o','d'};

figure('Position',[680 546 662 332])
hold on

if out_type == 1
    boxchart(((data_filtered.auto_corr_indx)),data_filtered.avg_benefit,'BoxWidth',1,'MarkerStyle','.','MarkerSize',7);    
    for i = 1:4
        plot(1:2:length(one_step_auto_corr),temp_benefit(:,i,1),'Marker',marker_type{i},'MarkerSize',10,'LineStyle','none','LineWidth',1.5);
    end
    ylabel('Avg fitness')
else
    boxchart(((data_filtered.auto_corr_indx)),data_filtered.var_benefit,'BoxWidth',1,'MarkerStyle','.','MarkerSize',7);
    for i = 1:4
        plot(1:2:length(one_step_auto_corr),temp_benefit(:,i,2),'Marker',marker_type{i},'MarkerSize',10,'LineStyle','none','LineWidth',1.5);
    end
    ylabel('Var fitness')
end

grid on
ax = gca;
ax.FontSize = 14;
ax.XTick = 1:2:length(one_step_auto_corr);
ax.XTickLabel = one_step_auto_corr(1:2:end);
xlim([0 38])
xlabel('\rho_{xx}(1)');
if out_type == 1
    saveas(gcf,['Avg_Benefit_vs_auto_corr_varied_m_and_c0_q_' num2str(q) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
else
    saveas(gcf,['Var_Benefit_vs_auto_corr_varied_m_and_c0_q_' num2str(q) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
end
end
%% similarity between probability of env and state transitions
optimal_c_N = zeros(length(one_step_auto_corr),2);
figure('Position',[183 538 907 340])
color_seq = {'m','cost','\rho_{xx}(1)'};
subplot(1,2,1)
for i = 1:length(one_step_auto_corr)
    % scatter(data_table{data_table.auto_corr_indx == one_step_auto_corr(i),9},data_table{data_table.auto_corr_indx == one_step_auto_corr(i),8},20,data_table{data_table.auto_corr_indx == one_step_auto_corr(i),3},'filled')
    scatter(data_table{data_table.auto_corr_indx == i,9},data_table{data_table.auto_corr_indx == i,8},20,one_step_auto_corr(data_table{data_table.auto_corr_indx == i,3}),'filled')

    hold on
    indices = find(data_table.auto_corr_indx == i);
    [~,max_indx] = max(data_table{indices,4});
    plot(data_table{indices(max_indx),9},data_table{indices(max_indx),8},'*','MarkerSize',10,'Color','r');
    optimal_c_N(i,:) = [data_table{indices(max_indx),1} c0_set(data_table{indices(max_indx),2})];

end

x = 0:0.1:1;
y = x;
plot(x,y,'--k','LineWidth',2);
grid on
ylabel('t_{S^l \rightarrow S^h}')
xlabel('t_{E^l \rightarrow E^h}')
ax= gca;
axis square
ax.FontSize = 14;
ax.XTick = 0:0.2:1;
ax.YTick = 0:0.2:1;
cb = colorbar;
cb.Ticks = one_step_auto_corr(1:6:end);

title(cb,'\rho_{xx}(1)')

subplot(1,2,2)
for i = 1:length(one_step_auto_corr)
    scatter(data_table{data_table.auto_corr_indx == i,11},data_table{data_table.auto_corr_indx == i,10},20,one_step_auto_corr(data_table{data_table.auto_corr_indx == (i),3}),'filled')
    hold on
    indices = find(data_table.auto_corr_indx == i);
    [~,max_indx] = max(data_table{indices,4});
    plot(data_table{indices(max_indx),11},data_table{indices(max_indx),10},'*','MarkerSize',10,'Color','r');
end

x = 0:0.1:1;
y = x;
plot(x,y,'--k','LineWidth',2);
grid on
axis square

ylabel('t_{S^h \rightarrow S^l}')
xlabel('t_{E^h \rightarrow E^l}')
ax= gca;
ax.FontSize = 14;
ax.XTick = 0:0.2:1;
ax.YTick = 0:0.2:1;
cb = colorbar;
cb.Ticks = one_step_auto_corr(1:6:end);

title(cb,'\rho_{xx}(1)')
saveas(gcf,['Transition_prob_comparison_varied_m_and_c0_q_' num2str(q) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])

figure('Position',[680 527 411 351])
tiledlayout(2,1,"TileSpacing","compact",Padding="compact")

nexttile(1)
plot(one_step_auto_corr,reshape(optimal_c_N(:,1),length(one_step_auto_corr),1),'LineWidth',2);
hold on
ylabel('Optimal m')
% xlabel('\rho_{xx}(1)')
grid on
ax =gca;
ax.FontSize = 14;
ax.XTick = one_step_auto_corr(1:6:end);
ax.YTick = N(1:2:end);
ylim([min(N) max(N)]);
nexttile(2)
plot(one_step_auto_corr,reshape(optimal_c_N(:,2),length(one_step_auto_corr),1),'LineWidth',2);
hold on
ylabel('Optimal cost')
xlabel('\rho_{xx}(1)')
ax =gca;
grid on
ax.FontSize = 14;
ax.XTick = one_step_auto_corr(1:6:end);
ax.YTick = [0 round(c0_set(10:10:end),2)];
ylim([0 max(c0_set)]);

saveas(gcf,['Optimal_m_and_c0_q_' num2str(q) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])
%% Variation of E[En Sn] with auto_corr
color_code = ["#D95319", "#EDB120", "#7E2F8E", 	"#77AC30", 	"#0072BD"];
indices = zeros(length(one_step_auto_corr),2);
for corr_indx = 1:length(one_step_auto_corr)
    indices(corr_indx,:) = [find(N == optimal_c_N(corr_indx,1)) find(c0_set == optimal_c_N(corr_indx,2))];
    optimal_corr(corr_indx) = env_cell_state_corr(indices(corr_indx,2),corr_indx,indices(corr_indx,1));
end

figure('Position',[680 543 365 335])
plt_indx = 1;
for c0_indx = [1 length(c0_set)]
    for N_indx = [1 length(N)]
        data3 = reshape(env_cell_state_corr(c0_indx,:,N_indx),1,length(one_step_auto_corr));
        plot(one_step_auto_corr,data3,LineWidth=2,Color=color_code(plt_indx));
        hold on
        plt_indx = plt_indx+1;
    end
end
plot(one_step_auto_corr,optimal_corr,'--',LineWidth=3,Color=color_code(end));
ylabel('Corr(E_n, S_n)')
xlabel('\rho_{xx}(1)')
% axis square
grid on
ax =gca;
ax.FontSize = 14;
ax.XTick = one_step_auto_corr(1:6:end);
ylim([-0.45 0.85])
saveas(gcf,['Cross_corr_env_cell_state_varied_m_and_c0_q_' num2str(q) '_p_' num2str(p) ' match reward ' num2str(reward) ' mismatch reward ' num2str(cost) '.png'])