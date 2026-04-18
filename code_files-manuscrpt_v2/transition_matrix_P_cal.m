% initialization of transition probability matrix
function P = transition_matrix_P_cal(N,p)
P = zeros(N+1,N+1);
for i = 1:N+1
    for j = 1:N+1
        if(i == 1 && j == 1)
            P(i,j) = 1-p;
        elseif(i == 1 && j == 2)
            P(i,j) = p;
        elseif(i> 1 && i <=N && (j-i == -1))
            P(i,j) =  (i-1)*(1-p)/N;
        elseif(i> 1 && i <=N && (j-i == 0))
            P(i,j) = (i-1)*p/N + (N-(i-1))*(1-p)/N;
        elseif(i> 1 && i <=N && (j-i == 1))
            P(i,j) = (N-(i-1))*p/N;
        elseif(i == N+1 && j == N)
            P(i,j) = 1-p;
        elseif(i == N+1 && j == N+1)
            P(i,j) = p;
        end

        % if(i == 1 && j == 1)
        %     P(i,j) = 1-p;
        % elseif(i == 1 && j == 2)
        %     P(i,j) = p;
        % elseif(i> 1 && i <=N && (j-i == -1))
        %     P(i,j) =  p*(1-p);
        % elseif(i> 1 && i <=N && (j-i == 0))
        %     P(i,j) = p^2+ (1-p)^2;
        % elseif(i> 1 && i <=N && (j-i == 1))
        %     P(i,j) = p*(1-p);
        % elseif(i == N+1 && j == N)
        %     P(i,j) = 1-p;
        % elseif(i == N+1 && j == N+1)
        %     P(i,j) = p;
        % end

%         if(i == 1 && j == 1)
%             P(i,j) = (1-p)*P_k_n_cal(N,i-1,p);
%         elseif(i == 1 && j == 2)
%             P(i,j) = p*P_k_n_cal(N,i-1,p);
%         elseif(i> 1 && i <=N && (j-i == -1))
%             P(i,j) =  P_k_n_cal(N,i-1,p)*(i-1)*(1-p)/N;
%         elseif(i> 1 && i <=N && (j-i == 0))
%             P(i,j) = P_k_n_cal(N,i-1,p)*((i-1)*p/N + (N-(i-1))*(1-p)/N);
%         elseif(i> 1 && i <=N && (j-i == 1))
%             P(i,j) = P_k_n_cal(N,i-1,p)*((N-(i-1))*p/N);
%         elseif(i == N+1 && j == N)
%             P(i,j) = P_k_n_cal(N,i-1,p)*(1-p);
%         elseif(i == N+1 && j == N+1)
%             P(i,j) = P_k_n_cal(N,i-1,p)*p;
%         end
    end
end
end

% function P_k_n = P_k_n_cal(N,k,p)
%     N
%     p
%     k
%     P_k_n = nchoosek(N+1,k)*(p^k)*(1-p)^(N+1-k);
% end