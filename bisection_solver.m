function x = bisection_solver(fun, x_left, x_right)
    if sign(fun(x_left)) == sign(fun(x_right))
        error('The initial interval does not cross zero. Bisection method will not work!');
    end
   
    max_it = 1000;
    a_thresh = 1e-14;
    b_thresh = 1e-14;

    c_old = x_left;

    for i = 1:max_it
        c = (x_left + x_right) / 2;
        
        delta_x = abs(c - c_old);
        f_c = fun(c);

        if delta_x < a_thresh || abs(f_c) < b_thresh
            x = c;
            fprintf('Converged after %d iterations.\n', i);
            return;
        end
        
        if sign(fun(x_left)) == sign(f_c) 
            x_right = c;
        else
            x_left = c;
        end
        
        c_old = c;
    end


end