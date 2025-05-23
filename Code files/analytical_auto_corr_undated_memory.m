
function auto_corr = analytical_auto_corr(N,p,max_lag)
P = transition_matrix_P_cal(N,p);
lags = 0:max_lag;

E_n_k_given_n = zeros(N+1,length(lags)); % E[X(n+k)|X(n) = a]
for k = lags
    for a = 1:N+1
        x_init = zeros(1,N+1);
        x_init(a) = 1;
        E_n_k_given_n(a,k+1) = x_init*(P^k)*([0:N]');
    end
end

x_eqb = x_init*P^(100);
E_n = x_eqb*[0:N]'; % E[X(n)]
E_n_2 = x_eqb*([0:N]').^2; % E[X(n)^2]
var_n = E_n_2 - (E_n)^2; % Var(X(n))

auto_corr = zeros(length(lags),1);

for k = lags
    auto_corr_sum = 0;
    for a = 1:N+1
        auto_corr_sum = auto_corr_sum + E_n_k_given_n(a,k+1)*x_eqb(a)*(a-1);
    end
    auto_corr(k+1) = (auto_corr_sum - E_n^2)/var_n;
end
end
% plot(0:N+1,auto_corr);