function env = correlated_env(p1,p2,total_time)

env = zeros(1,total_time);
env(1) = binornd(1,(p1/(p1+p2)));
% if(~isempty(varargin))
%     env(1) = varargin(1);
% else
%     if(p2 == 0)
%         env(1) = 1;
%     elseif(p1 == 0)
%         env(1) = 0;
%     else
%         env(1) = randi(2)-1;
%     end
% end
T = [1-p1;p1; p2; 1-p2];

rand_num = rand(total_time,1);

for t = 1:total_time-1
    a = rand_num(t);
    if(env(t) == 0)
        for j = 1:2
            if(sum(T(1:j)) > a)
                rxn = j;
                break;
            end
        end
    else
        for j = 3:4
            if(sum(T(3:j)) > a)
                rxn = j;
                break;
            end
        end
    end

    switch rxn
        case {1,3}
            env(t+1) = 0;
        case {2,4}
            env(t+1) = 1;
    end
end
end