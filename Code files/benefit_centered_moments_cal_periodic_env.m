%%%% calculation of centered momenets of cellular benefit in periodic
%%%% environments

function [centered_moments_benefit, centered_moments_benefit_with_time] = benefit_centered_moments_cal_periodic_env(T,r,c, prob_state_1, prob_state_1_given_1,cost)

r1 = r(1);
r0 = r(2);
c0 = c(2);
c1 = c(1);
centered_moments_benefit = zeros(1,1);
uncentered_moments_benefit = zeros(T,2);

for m = 1:T
    if(m <= floor(T/2))
        Em = 0;
    else
        Em = 1;
    end

    if(m > 2)
%     uncentered_moments_benefit(m,1) = prob_state_1(m-1)*(Em*r1+(1-Em)*c1) + (1-prob_state_1(m-1))*((1-Em)*r0+Em*c0) + 2*cost*(prob_state_1_given_1(m-1) - prob_state_1(m-1));
    uncentered_moments_benefit(m,1) = prob_state_1(m-1)*(Em*r1+(1-Em)*c1) + (1-prob_state_1(m-1))*((1-Em)*r0+Em*c0) + cost*(prob_state_1_given_1(m-1)*prob_state_1(m-2) - prob_state_1(m-1)) + cost*(prob_state_1_given_1(m-1)*prob_state_1(m-2) - prob_state_1(m-2));
%     uncentered_moments_benefit(m,2) = (1-prob_state_1(m-1))*((1-Em)*r0^2+Em*c0^2) + prob_state_1(m-1)*(Em*r1^2+(1-Em)*c1^2) + 2*(prob_state_1(m-1) - prob_state_1_given_1(m-1))*(cost^2  - cost*((1-Em)*r0+Em*c0) - cost*(Em*r1+(1-Em)*c1));
    elseif(m == 1)
    uncentered_moments_benefit(m,1) =  prob_state_1(end)*(Em*r1+(1-Em)*c1) + (1-prob_state_1(end))*((1-Em)*r0+Em*c0) + cost*(prob_state_1_given_1(end)*prob_state_1(end-1) - prob_state_1(end)) + cost*(prob_state_1_given_1(end)*prob_state_1(end-1) - prob_state_1(end-1));
%     uncentered_moments_benefit(m,2) =  (1-prob_state_1(end))*((1-Em)*r0^2+Em*c0^2) + prob_state_1(end)*(Em*r1^2+(1-Em)*c1^2) + 2*(prob_state_1(end) - prob_state_1_given_1(end))*(cost^2  - cost*((1-Em)*r0+Em*c0) - cost*(Em*r1+(1-Em)*c1));
    elseif(m == 2)
        uncentered_moments_benefit(m,1) =  prob_state_1(1)*(Em*r1+(1-Em)*c1) + (1-prob_state_1(1))*((1-Em)*r0+Em*c0) + cost*(prob_state_1_given_1(1)*prob_state_1(end) - prob_state_1(1)) + cost*(prob_state_1_given_1(1)*prob_state_1(end) - prob_state_1(end));
    end
end

centered_moments_benefit_with_time = uncentered_moments_benefit(:,1);
centered_moments_benefit(1) = sum(uncentered_moments_benefit(:,1))/T;
% centered_moments_benefit(2) = sum((uncentered_moments_benefit(:,2) - (uncentered_moments_benefit(:,1)).^2))/T^2;

end