% initialization of State Transition Matrix for undated memory

function P = transition_matrix_P_cal_corr_env(N,p1,p2)
P = zeros(2*(N+1),2*(N+1));

kn_0 = 0:N;
kn_1 = 0:N;
kn = [kn_0 kn_1];
for i = 1:2*(N+1)
    for j = 1:2*(N+1)

        if(i < N+1)
            trans_prob_vec = [(kn(i)/N)*(1-p1) (1-kn(i)/N)*(1-p1)  (kn(i)/N)*(p1) (1-kn(i)/N)*(p1)];
            if(j == i-1)
                P(i,j) = trans_prob_vec(1);
            elseif(j == i)
                P(i,j) = trans_prob_vec(2);
            elseif(j == i+(N+1))
                P(i,j) = trans_prob_vec(3);
            elseif(j == i+(N+1)+1)
                P(i,j) = trans_prob_vec(4);
            end
        elseif(i == N+1)
            P(N+1,N+1) = 1;
        elseif(i == N+2)
            P(N+2,N+2) = 1;
        else
            trans_prob_vec = [(kn(i)/N)*(p2) (1-kn(i)/N)*(p2)  (kn(i)/N)*(1-p2) (1-kn(i)/N)*(1-p2)];
            if(j == i-(N+1)-1)
                P(i,j) = trans_prob_vec(1);
            elseif(j == i-(N+1))
                P(i,j) = trans_prob_vec(2);
            elseif(j == i)
                P(i,j) = trans_prob_vec(3);
            elseif(j == i+1)
                P(i,j) = trans_prob_vec(4);
            end
        end
    end
end

% restructing the transition matrix 
P_temp = zeros(2*(N+1),2*(N+1));
P_temp1 = zeros(2*(N+1),2*(N+1));

for i = 1:N+1
       P_temp(2*(i-1)+1,:) = P(i,:);
       P_temp(2*(i-1)+2,:) = P(i+N+1,:);
end

for j = 1:N+1
       P_temp1(:,2*(j-1)+1) = P_temp(:,j);
       P_temp1(:,2*(j-1)+2) = P_temp(:,j+N+1);
end

P = P_temp1;

end