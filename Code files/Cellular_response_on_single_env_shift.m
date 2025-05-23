close all
clear all
results = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\Phenotypic_adaptation_Costly_transition\Figures for manuscript';
code_files = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\Phenotypic_adaptation_Costly_transition\Code files';
% models = {'Homo MC', 'Inhomo MC', 'Original'};

q = 0.3; % indifference environment 

alpha = 1; % difference between (r_h^+)-(r_l^-) 
beta = alpha*q./(1-q); % (r_l^+)-(r_h^-)
cost = [1 1]; % (r_h^-), (r_l^-)
reward = [alpha+cost(2) beta+cost(1)]; % (r_h^+), (r_l^+)

c0 = linspace(0,min(alpha,beta)-min(alpha,beta)/100,3);


line_style = {'-','--'};
color_order = ["#0072BD", "#D95319",  "#77AC30", "#EDB120", "#7E2F8E"];
% color_order = ["#0072BD",  "#77AC30", "#7E2F8E"];

% benefit of mismatching the env is considered 1

N = 8;%2:2:10; % memory size

p2 = {[q-0.1];
    [q+0.1]};

markers = ['o','x','d','^','+','>'];
num_time_steps = 200;

for q_indx = 2
    if(q_indx == 1)
        p1 = 1;
    else
        p1 =0;
    end

    if(~isempty(p2{q_indx}))
        adapt_time = zeros(length(N),length(c0),2);
        adapt_prob_with_time = zeros(length(N),length(c0),2,num_time_steps);
        avg_benefit_with_time = zeros(length(N),length(c0),2,num_time_steps);
        prob_resi_optimal_state = zeros(length(N),length(c0),2,num_time_steps);

        for c0_indx = 1:length(c0)

            q01 = (beta+c0(c0_indx))/(alpha+beta);
            q10 = (beta-c0(c0_indx))/(alpha+beta);
            q_cost = [q10 q01];

            for N_indx = 1:length(N)
                for model_indx = 1:2
                    [adapt_time(N_indx,c0_indx,model_indx), prob_resi_optimal_state(N_indx,c0_indx,model_indx,:),avg_benefit_with_time(N_indx,c0_indx,model_indx,:),adapt_prob_with_time(N_indx,c0_indx,model_indx,:)] = adaptation_func(model_indx,[p1 p2{q_indx}(1)],q_cost(q_indx),N(N_indx),alpha,beta,num_time_steps,q_cost);
                end
            end
        end
    end
end

%% Cell state and average benefit time trajectories with cost of transition variation
N_indx = 1;

color_order = ["#0072BD", "#77AC30", "#D95319"];

figure('Position',[301.8000 317.8000 945.6000 444.2000]);
% subplot(2,5,N_indx)
tiledlayout(2,2);
for model_indx = 1:2
    nexttile
    plot(1:num_time_steps,reshape(prob_resi_optimal_state(N_indx,:,model_indx,:),length(c0),num_time_steps),'LineWidth',2,'LineStyle',line_style{model_indx});

    if(N_indx == 1)
        ylabel('P(S_n = 1)')
        lg = legend(num2str(round(c0,2)'),'Location','southeast');
        title(lg,'c_0')
    end

    ax = gca;
    xlabel('time steps')
    ax.FontSize = 16;
    xlim([0 40])
    ax.XTick = 0:10:num_time_steps;
    grid on
    ylim([0 1]);
    ax.YTick = 0:0.2:1;

    colororder(color_order)
end

% saveas(gcf, ['Prob residing in optimal state c0 variation  q ' num2str(q) ' p2 ' num2str(p2{q_indx}) '.png']);

figure('Position',[301.8000 317.8000 945.6000 444.2000]);
% subplot(2,5,N_indx)
tiledlayout(2,2);
for model_indx = 1:2
nexttile
plot(1:num_time_steps,reshape(avg_benefit_with_time(N_indx,:,model_indx,:),length(c0),num_time_steps),'LineWidth',2,'LineStyle',line_style{model_indx});
hold on
colororder(color_order)
plot(1:num_time_steps,reshape(avg_benefit_with_time(N_indx,:,2,:),length(c0),num_time_steps),'LineWidth',2,'LineStyle',line_style{2});

% axis square
% title({['q ' num2str(q) ', p_2 ' num2str(p2{q_indx})]}, {[ 'p_1 ' num2str(p1) ', N ' num2str(N(N_indx))]})

if(N_indx == 1)
    ylabel('Average Benefit')
    lg = legend(num2str(round(c0,2)'),'Location','southeast');
    title(lg,'c_0')
end
ax = gca;
xlabel('time steps')
xlim([0 num_time_steps])
ylim([min(avg_benefit_with_time(N_indx,:,:,:),[],"all")-0.01 max(avg_benefit_with_time(N_indx,:,:,:),[],"all")+0.01])
ax.XTick = 0:10:num_time_steps;
ax.FontSize = 16;
grid on
end

% saveas(gcf, ['Average benefit with time c0 variation q ' num2str(q) ' p2 ' num2str(p2{q_indx}) '.png']);

%% Mean 0 -> 1 switch time for increasing cost

models = {'Undated', 'Dated'};

alpha = 1;% alpha = r1 - c0
q = 0.3;%0.1:0.1:0.9;6
beta = alpha*q./(1-q); % beta = r0 - c1

c0 = linspace(0,min(alpha,beta)-min(alpha,beta)/100,500);


line_style = {'-','--'};
color_order = ["#0072BD", "#D95319",  "#77AC30","#7E2F8E"];
% benefit of mismatching the env is considered 1

N = 8;
% N = 5;
% p2 = 0:0.1:q-0.1;
% p2 = {0:0.05:q-0.05;
%     1:-0.05:q+0.05};

p2 = [0.4 0.7 1];
p1 = 0;%:0.01:1;

num_time_steps = 50;

centered_moments = zeros(length(p2),length(c0),2,2);
adapt_time = zeros(length(p2),length(c0),2);
state_difference = zeros(length(p2),length(c0),2);


for p2_indx = 1:length(p2)
    parfor c0_indx = 1:length(c0)
        q01 = (beta+c0(c0_indx))/(alpha+beta);
        q10 = (beta-c0(c0_indx))/(alpha+beta);
        q_cost = [q10 q01];
        [temp_centered_moments_1,~,~,~,~] = env_switch_mean_adaptation_time_undated_memory(N,q_cost(2),q_cost(1),p1,p2(p2_indx));
        [temp_centered_moments_2,~,~,~,~] = env_switch_mean_adaptation_time_dated_memory(N,q_cost(2),q_cost(1),p1,p2(p2_indx));
        centered_moments(p2_indx,c0_indx,:,:) = [temp_centered_moments_1(1,:)' temp_centered_moments_2(1,:)'];

%         [adapt_time_1,state_difference_1] = env_switch_undated_memory_cell_state_eqb_time_and_state_diff(N,q_cost(2),q_cost(1),p1,p2(p2_indx));
%         [adapt_time_2,state_difference_2] = env_switch_dated_memory_cell_state_eqb_time_and_state_diff(N,q_cost(2),q_cost(1),p1,p2(p2_indx));

%         adapt_time(p2_indx,c0_indx,:) = [adapt_time_1 adapt_time_2];
%         state_difference(p2_indx,c0_indx,:) = [state_difference_1 state_difference_2];
    end
end


figure('Units','normalized','Position',[0 0 1 1])
subplot(2,3,1)

line_type = {'-','--'};
color_order = ["#0072BD","#77AC30", "#D95319"];
for p2_indx = 1:length(p2)
    for model_indx = 1:2
        % Plot the errorbars
        plot(c0,centered_moments(p2_indx,:,1,model_indx)',line_type{model_indx},'LineWidth',2,'Color',color_order(p2_indx));
        hold on
    end
end

hold on
set(gca,'YScale','log')
ax = gca;
% ax.XTickLabel = c0_set;
ax.FontSize = 16;
grid on
xlabel('Transition cost (c_{lh} = c_{hl})')
ylabel(['Mean S^l to S^h switch time']);
% legend(['\tau_' state_label{1} '^A'], ['\tau_' state_label{1} '^N'],['\tau_' state_label{2} '^A'], ['\tau_' state_label{2} '^N'])
% colororder(color_order)
xlim([min(c0) max(c0)])

% legend('Undated_A','Undated_N','Dated_A','Dated_N')
saveas(gcf,['Mean S^l to S^h switch time q ' num2str(q) ' N ' num2str(N) ' p1 ' num2str(p1) '.png']);


function [adapt_time, prob_resi_optimal_state,avg_benefit_with_time,adapt_prob_with_time] = adaptation_func(model,p,q,N,alpha,beta,num_time_steps,q_cost)

if(model == 1)
    [centered_moments,prob_resi_time_state_0,prob_resi_time_state_1,prob_state_1_start_0,prob_state_1_start_1] = env_switch_mean_adaptation_time_undated_memory(N,q_cost(2),q_cost(1),p(1),p(2),num_time_steps);
else
    [centered_moments,prob_resi_time_state_0,prob_resi_time_state_1,prob_state_1_start_0,prob_state_1_start_1] = env_switch_mean_adaptation_time_dated_memory(N,q_cost(2),q_cost(1),p(1),p(2),num_time_steps);
end
if(p(1) == 0)
    adapt_time = centered_moments(1,1);
else
    adapt_time = centered_moments(1,2);
end

reward = [alpha+1 beta+1];
cost= [1 1];
S_1 =p(2)*reward(1) + (1-p(2))*cost(1);
S_0 =(1-p(2))*reward(2) + (p(2))*cost(2);

n = 0;
adapt_prob_with_time = zeros(num_time_steps,1);
avg_benefit_with_time = zeros(num_time_steps,1);
prob_resi_optimal_state = zeros(num_time_steps,1);
while(n < num_time_steps)
    n = n+1;
    if(p(1) == 0)
        adapt_prob_with_time(n) = sum(prob_resi_time_state_0(1:n));
        prob_resi_optimal_state(n) = prob_state_1_start_0(n);
        avg_benefit_with_time(n) = prob_resi_optimal_state(n)*S_1+(1-prob_resi_optimal_state(n))*S_0;

    else
        adapt_prob_with_time(n) = sum(prob_resi_time_state_1(1:n));
        prob_resi_optimal_state(n) = prob_state_1_start_1(n);
        avg_benefit_with_time(n) = (1-prob_resi_optimal_state(n))*S_1+(prob_resi_optimal_state(n))*S_0;
    end
end

end