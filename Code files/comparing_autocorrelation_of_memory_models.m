%% Autocorrelation function
% color_order = ["#0072BD", "#D95319", "#EDB120","#7E2F8E"];
color_order = ["#000000","#ADA1A1"];

% fig1 = figure('Units','normalized','Position',[0 0 1 1]);
% for plotting autocorrelation of cell memory and state for the two models across p values
line_type = {'-','--'};
N = [4:4:12];
p= 0.4;
figure('Position',[680 646 339 332]);
tiledlayout(length(N),1,TileSpacing="tight")
for N_indx = 1:length(N)
    nexttile
    for model_indx = 1:2
        % undated memory model
        if(model_indx == 1)
            lags = N(end)+1;
            autocorr_memory = analytical_auto_corr_undated_memory(N(N_indx),p,lags);
        else
            autocorr_memory = [1-(0:N(N_indx))/N(N_indx) zeros(1,length(autocorr_memory)-(N(N_indx)+1))];
        end
        %         plot(0:length(autocorr_memory)-1, autocorr_memory,'Marker',marker_type{model_indx},'MarkerSize',5,'LineWidth',2,Color=color_order(N_indx));
        plot(0:length(autocorr_memory)-1, autocorr_memory,'LineStyle',line_type{model_indx},'LineWidth',2,Color=color_order(1));
        hold on
    end
    grid on
    % ylabel('Autocorrelation (cell memory)')
    ax = gca;
    ax.FontSize = 11;
    ax.YTick = 0:0.2:1;
    ylim([0 1])
    xlim([0 N(end)+2]);
end
xlabel('Lags')

saveas(gcf,['Numerical Autocorrelation of memory models N ' num2str(N) ' p ' num2str(p) '.png']);