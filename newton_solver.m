

%Definition of the test function and its derivative
% a_thresh = 1e-14;
% b_thresh = 1e-14;


% function root = newton_method(fun, a, b, a_thresh, b_thrsh, max_iter)
%     delta_x = abs()
%     while (delta_x > a_thresh && )
% end

% test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
% test_derivative01 = @(x) 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
% 
% 
% [root, guesses] = newton_solver(test_func01, test_derivative01, 091, a_thresh, b_thresh, 200); 


function [root, guesses] = newton_solver(fun, dfun, x0) 

    max_iter = 1000;
    x_thresh = 1e-14;
    y_thresh = 1e-14;

    i = 0;
    guesses = [x0];
    while (i < max_iter)
        x_n = x0 - fun(x0)/dfun(x0);
        guesses(end + 1) = x_n;
        root = x_n;
        if (abs(x_n - x0) > x_thresh || abs(fun(x0)) > y_thresh)
            i = i + 1;
        else 
            i = max_iter + 1;
        end
        x0 = x_n;
    end
    guesses(end + 1) = root;
end
