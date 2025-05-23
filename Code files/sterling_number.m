function out = sterling_number(n,k)
% sterling number of second kind
sum = zeros(n+1,1);
for j = 0:n
    for i = 0:j
        sum(j+1) = sum(j+1) + ((-1)^(j-i)*i^n)/(factorial(j-i)*factorial(i));
    end
end
out = sum(k+1);
end