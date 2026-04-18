function env = discrete_period_biased_env(T,total_duration,Eh_bias)

El_bias = 1-Eh_bias;
T_El = floor(T*El_bias);
T_Eh = T - T_El;

env = zeros(total_duration,1);
cycle = 1;

for indx = 1:total_duration
    if(indx <= T_Eh + (cycle-1)*T)
        env(indx) = 1;
    else
        env(indx) = 0;
    end

    if(indx == cycle*T)
        cycle = cycle + 1;
    end
end

end