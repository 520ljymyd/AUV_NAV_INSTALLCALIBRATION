function error = optimization_function(params)
n_1 = round(params(1)/5)*5;  % n_1 is the first parameter
n_2 = round(params(2)/5)*5;  % n_2 is the second parameter


% Ensure that n_1 > n_2, if not, swap them
if n_1 <= n_2
    temp = n_1;
    n_1 = n_2;
    n_2 = temp;
end

angle = func_angle(n_1,n_2);
error = abs(mean(angle)-0.2);
end