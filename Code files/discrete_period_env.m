function env = discrete_period_env(T,total_duration)

env = zeros(total_duration,1);
cycle = 1;

for indx = 1:total_duration
    if(indx <= (T/2)+ (cycle-1)*T)
        env(indx) = 0;
    else
        env(indx) = 1;
    end

    if(indx == cycle*T)
        cycle = cycle + 1;
    end
end

end