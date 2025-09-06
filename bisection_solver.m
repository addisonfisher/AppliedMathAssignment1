function [x, guesses_c, guesses_fc, guesses_it] = bisection_solver(fun, x_left, x_right)
    if sign(fun(x_left)) == sign(fun(x_right))
        error('The initial interval does not cross zero. Bisection method will not work!');
    end
   
    max_it = 1000;
    a_thresh = 1e-14;
    b_thresh = 1e-14;

    c_old = x_left;
    
    guesses_c = [];
    guesses_fc = [];
    guesses_it = [];

    for i = 1:max_it
        c = (x_left + x_right) / 2;
        
        delta_x = abs(c - c_old);
        f_c = fun(c);

        if i > 0
            guesses_c(end+1) = c_old;
            guesses_fc(end+1) = fun(c_old);
            guesses_it(end+1) = i;
        end

        if delta_x < a_thresh || abs(f_c) < b_thresh
            x = c;
            fprintf('Converged after %d iterations.\n', i);
            return;
        end
        
        if sign(fun(x_left)) == sign(f_c) 
            x_left = c;
        else
            x_right = c;
        end
        
        c_old = c;
    end


end