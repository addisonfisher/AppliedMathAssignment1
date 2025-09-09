%Secant Solver Problem

function [xn_list, xprev_list, x_count, x_root] = secant_method(fun, x0, x1, a_thresh, b_thresh, max_iter)

    %Secant solver takes two initial inputs and updates new guess off of
    %the secant line of the two inputs
    
    xn_list = [];
    xprev_list = [];
    x_count = [];

    x_n = x1;
    x_prev = x0; 
    
    count = 1; 

    for i = max_iter
        
        %Defining function at guessed root
        [f_n, ~] = fun(x_n);
        
        %Defining function at previous guessed root
        [f_prev, ~] = fun(x_prev);
        
        if abs(f_n - f_prev) == 0
            disp ("Update step is too large")
            break
        end 
        %Update new root for secant method
        x_next = x_n - f_n * ((x_n - x_prev)/(f_n - f_prev));
        
        if abs(x_next - x_prev) < a_thresh && f_n < b_thresh
            disp("Root found")
            break
        end

        %Update previous root for secant method 
        xprev_list = [xprev_list, x_prev];
        x_prev = x_n; 
        
        %Update new root for secant method
        x_n = x_next;
        xn_list = [xn_list, x_n];

        x_count = [x_count, count];
        count = count + 1; 
    end 

    if i == max_iter
        disp("Root not found")
    end

    x_root = x_n;
end 