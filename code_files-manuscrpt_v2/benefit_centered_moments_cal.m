%%%% calculation of centered momenets of cellular benefit in stochastic
%%%% environments

function [centered_moments_benefit,env_cell_state_corr, P_Sn_En] = benefit_centered_moments_cal(p1,p2,r, c, cond_prob_Sn_Sn_1, stat_dist_Sh, stat_dist_Sl,stat_dist_Sn_Sn_1_Sh,cost)

rh = r(1);
rl = r(2);
cl = c(2);
ch = c(1);
centered_moments_benefit = zeros(2,1);

% stationary distribution of memory states
prob_sum = zeros(1,4);
% prob_sum(1) =  P[Sn = Sh, En-1 = Eh]
% prob_sum(2) =  P[Sn = Sh, En-1 = El]
% prob_sum(3) =  P[Sn = Sh, En-1 = Eh | Sn-1 = Sh]
% prob_sum(4) = P[Sn = Sh, En-1 = El | Sn-1 = Sh]

for i = 1:length(stat_dist_Sh)
    if(mod(i,2) == 0)
        prob_sum(1) = prob_sum(1) + stat_dist_Sh(i);
        prob_sum(3) = prob_sum(3) + stat_dist_Sn_Sn_1_Sh(i);
    else
        prob_sum(2) = prob_sum(2) + stat_dist_Sh(i); % P[Sn = Sh, En-1 = El]
        prob_sum(4) = prob_sum(4) + stat_dist_Sn_Sn_1_Sh(i);
    end
end


E_En = p1/(p1+p2);
E_Sn = sum(stat_dist_Sh);
E_En_Sn = (1-p2)*prob_sum(1) + (p1)*prob_sum(2);
E_Sn_Sn_1 = cond_prob_Sn_Sn_1*E_Sn;
E_En_Sn_1 = (prob_sum(1)*((1-p2)^2+p1*p2)+prob_sum(2)*(p1*(1-p1)+p1*(1-p2)));
E_En_Sn_Sn_1 = (prob_sum(3)*(1-p2)+prob_sum(4)*p1)*E_Sn;

centered_moments_benefit(1) = rl - E_Sn*(rl-ch) - (rl-cl)*E_En + (rl-cl+rh-ch)* E_En_Sn + 2*cost*(E_Sn_Sn_1 - E_Sn);

const1 = (cl^2-rl^2);
const2 = ch^2+ 2*cost^2 -rl^2-2*cost*(ch+rl);
const3 = -ch^2-cl^2+rh^2+rl^2+2*cost*(ch-rh);
const4 = -2*cost^2 + 2*cost*(ch+rl);
const5 = -2*cost*(cl-rl);
const6 = 2*cost*(-ch+cl+rh-rl); 

centered_moments_benefit(2) = rl^2 + E_En*const1 + E_Sn*const2 + E_En_Sn*(const3) + E_Sn_Sn_1*(const4) + E_En_Sn_1*(const5) + E_En_Sn_Sn_1*(const6) - centered_moments_benefit(1)^2;

std_Sn = sqrt(E_Sn*(1-E_Sn));
std_En = sqrt(E_En*(1-E_En));
env_cell_state_corr = (E_En_Sn - E_En*E_Sn)/(std_Sn*std_En);

%% Calculation of joint probability of cell and env state
P_Sn_Sh_given_En_1_Eh = (sum(stat_dist_Sh(2:2:end))/(sum(stat_dist_Sh(2:2:end))+sum(stat_dist_Sl(2:2:end))));
P_Sn_Sh_given_En_1_El = (sum(stat_dist_Sh(1:2:end))/(sum(stat_dist_Sh(1:2:end))+sum(stat_dist_Sl(1:2:end))));
P_Sn_Sl_given_En_1_Eh = (sum(stat_dist_Sl(2:2:end))/(sum(stat_dist_Sh(2:2:end))+sum(stat_dist_Sl(2:2:end))));
P_Sn_Sl_given_En_1_El = (sum(stat_dist_Sl(1:2:end))/(sum(stat_dist_Sh(1:2:end))+sum(stat_dist_Sl(1:2:end))));

P_Sh_Eh = P_Sn_Sh_given_En_1_Eh * (1-p2)* E_En + P_Sn_Sh_given_En_1_El *p1* (1-E_En);
P_Sh_El = P_Sn_Sh_given_En_1_Eh * (p2)* E_En + P_Sn_Sh_given_En_1_El *(1-p1)* (1-E_En);
P_Sl_Eh = P_Sn_Sl_given_En_1_Eh * (1-p2)* E_En + P_Sn_Sl_given_En_1_El *p1* (1-E_En);
P_Sl_El = P_Sn_Sl_given_En_1_Eh * (p2)* E_En + P_Sn_Sl_given_En_1_El *(1-p1)* (1-E_En);

P_Sn_En = [P_Sh_Eh P_Sh_El P_Sl_Eh P_Sl_El];

end