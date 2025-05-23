%%%% calculation of centered momenets of cellular benefit in stochastic
%%%% environments

function [centered_moments_benefit] = benefit_centered_moments_cal(p,r, c, centered_moments_resi_time, cond_prob_Sn_Sn_1,cost)

r1 = r(1);
r0 = r(2);
c0 = c(2);
c1 = c(1);
centered_moments_benefit = zeros(2,1);

%E(1-Sn)
mean_resi_state_0 = centered_moments_resi_time(1,1)/(sum(centered_moments_resi_time(1,:)));
%E(Sn)
mean_resi_state_1 = 1-mean_resi_state_0;

centered_moments_benefit(1) = mean_resi_state_1*(p*r1+(1-p)*c1) + mean_resi_state_0*((1-p)*r0+p*c0) + 2*cost*(cond_prob_Sn_Sn_1*mean_resi_state_1 - mean_resi_state_1);

% E(Sn^2) = E(Sn)
mom0 = mean_resi_state_0; %uncentered_moments_resi_time(2,1)/(uncentered_moments_resi_time(2,1)+uncentered_moments_resi_time(2,2));
mom1 = mean_resi_state_1; %1-mom0;
centered_moments_benefit(2) = mom0*((1-p)*r0^2+p*c0^2) + mom1*(p*r1^2+(1-p)*c1^2) + 2*(mom1 - cond_prob_Sn_Sn_1*mom1)*(cost^2  - cost*((1-p)*r0+p*c0) - cost*(p*r1+(1-p)*c1)) - centered_moments_benefit(1)^2;

% % Calculation of skewness
% %E[A^3]
% temp1= mean_resi_state_1*(p*r1^3+(1-p)*c1^3) + mean_resi_state_0 * ((1-p)*r0^3+p*c0^3);
% temp2 = -3*(mom1*(p*r1^2+(1-p)*c1^2) + mom0*((1-p)*r0^2+p*c0^2))*centered_moments_benefit(1);
% temp3 = 2*centered_moments_benefit(1)^3;
% 
% centered_moments_benefit(3) = (temp1+temp2+temp3)/(centered_moments_benefit(2)*sqrt(centered_moments_benefit(2)));

% temp_moments = zeros(4,2);
% for n = 1:4
%     for x = 1:length(prob_residence_time_state_0)
%         for y = 1:length(prob_residence_time_state_1)
%             temp_moments(n,1) = temp_moments(n,1) + ((x/(x+y))^n)*joint_prob_residence_time_state_0_and_1(x,y);
%             temp_moments(n,2) = temp_moments(n,2) + ((y/(x+y))^n)*joint_prob_residence_time_state_0_and_1(x,y);
%         end
%     end
% end
% 
% centered_momnets_resi_frac = [temp_moments(1,:); temp_moments(2,:)-temp_moments(1,:).^2];
% 
% % defining moments of benefit given a state for calculations of the moments
% 
% r0 = [beta+1,1]; % benefit on matching and mismatching of cell state 0 with the environment
% r1 = [alpha+1,1]; % benefit on matching and mismatching of cell state 1 with the environment
% moment_vec= zeros(14,1);
% p_vector_1 = [(1-p)*p; p^2; (1-p)^2; (1-p)*p];
% p_vector_2 = [p;0;(1-p);0];
% for mom_indx = 1:length(moment_vec)
%     temp_mat= zeros(2,2);
%     for i = 1:2
%         for j = 1:2
%             switch mom_indx
%                 case 1 %E[(a-b)]; where, RV 'a' take r0 and c0 values with prob 1-p and p, and RV 'b' take r1 (r1(1) and c1 (r1(2)) values with prob p and 1-p
%                     temp_mat(i,j) = (r0(i)-r1(j));
%                 case 2 %E[b]
%                     if(i == 1)
%                         temp_mat(i,j) = r1(j);
%                     end
%                 case 3 %E[(a-b)^2]
%                     temp_mat(i,j) = (r0(i)-r1(j))^2;
%                 case 4 %E[(a-b)b]
%                     temp_mat(i,j) = (r0(i)-r1(j))*r1(j);
% 
%                 case 5 %E[b^2]
%                     if(i == 1)
%                         temp_mat(i,j) = r1(j)^2;
%                     end
%                 case 6 %E[(a-b)^3]
%                     temp_mat(i,j) = (r0(i)-r1(j))^3;
%                 case 7 %E[(a-b)^2 b]
%                     temp_mat(i,j) = ((r0(i)-r1(j))^2)*r1(j);
%                 case 8 %E[(a-b) b^2]
%                     temp_mat(i,j) = (r0(i)-r1(j))*r1(j)^2;
%                 case 9 %E[b^3]
%                     if(i == 1)
%                         temp_mat(i,j) = r1(j)^3;
%                     end
%                 case 10 %E[(a-b)^4]
%                     temp_mat(i,j) = (r0(i)-r1(j))^4;
%                 case 11 %E[(a-b)^3 b]
%                     temp_mat(i,j) = ((r0(i)-r1(j))^3)*r1(j);
%                 case 12 %E[(a-b)^2 b^2]
%                     temp_mat(i,j) = ((r0(i)-r1(j))^2)*r1(j)^2;
%                 case 13 %E[(a-b) b^3]
%                     temp_mat(i,j) = ((r0(i)-r1(j))^1)*r1(j)^3;
%                 case 14 %E[b^4]
%                     if(i == 1)
%                         temp_mat(i,j) = r1(j)^4;
%                     end
%             end
%         end
%     end
% 
%     switch mom_indx
%         case {2,5,9,14}
%             moment_vec(mom_indx) = reshape(temp_mat,1,[])*p_vector_2;
%         otherwise
%             moment_vec(mom_indx) = reshape(temp_mat,1,[])*p_vector_1;
%     end
% end
% 
% mom = zeros(4,1);
% mom(1) = moment_vec(1)*temp_moments(1,1) + moment_vec(2);
% mom(2) = moment_vec(3)*temp_moments(2,1)+2*moment_vec(4)*temp_moments(1,1) + moment_vec(5);
% mom(3) = moment_vec(6)*temp_moments(3,1)+3*moment_vec(7) *temp_moments(2,1) + 3*moment_vec(8)*temp_moments(1,1) + moment_vec(9);
% mom(4) = moment_vec(10)*temp_moments(4,1)+4*moment_vec(11)*temp_moments(3,1) + 12*moment_vec(12)*temp_moments(2,1) + 4*moment_vec(13)*temp_moments(1,1) + moment_vec(14);
% 
% % calculation of mean (mu) E(X)
% centered_moments_benefit(1) = mom(1);
% % calculation of variance (sigma^2) E((X-mu)^2)
% centered_moments_benefit(2) = mom(2)-(mom(1))^2;
% % calculation of skewness E((X-mu)^3)/(sigma^3)
% centered_moments_benefit(3)= (mom(3) - 3*mom(2)*mom(1) + 2*mom(1)^3)/(centered_moments_benefit(2)*sqrt(centered_moments_benefit(2)));
% % calculation of kurtosis E((X-mu)^4)/(sigma^4)
% centered_moments_benefit(4)= (mom(4) - 4*mom(3)*mom(1) + 12*mom(2)*mom(1)^2 - 3*mom(1)^3)/(centered_moments_benefit(2)^2);
end