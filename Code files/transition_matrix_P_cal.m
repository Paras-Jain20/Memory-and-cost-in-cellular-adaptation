% initialization of State Transition Matrix for undated memory

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
    end
end
end