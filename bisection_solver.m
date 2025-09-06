function [x, guesses_x0, guesses_x1, guesses_it] = bisection_solver(fun, x_left, x_right)
    if sign(fun(x_left)) == sign(fun(x_right))
        error('The initial interval does not cross zero. Bisection method will not work!');
    end
   
    max_it = 1000;
    a_thresh = 1e-14;
    b_thresh = 1e-14;

    c_old = x_left;
    
    guesses_x0 = [];
    guesses_x1 = [];
    guesses_it = [];

    for i = 1:max_it
        c = (x_left + x_right) / 2;
        
        delta_x = abs(c - c_old);
        f_c = fun(c);

        if i > 0
            guesses_x0(end+1) = c_old;
            guesses_x1(end+1) = c;
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