function test

    solver_flag = 2;
    
    fun = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) - 0.7 - exp(x/6);
    
    % initial guess
    x_guess0 = -5;
    
    % guess lists
    num_trials = 250;
    guess_list1 = randi([-30, -5], 1, num_trials); % x_left values
    guess_list2 = randi([5, 30], 1, num_trials);    % x_right values
    
    % filtering constants
    filter_list = [1e-15, 1e-2, 1e-14, 1e-2, 2];
    
    convergence_analysis(solver_flag, fun, x_guess0, guess_list1, guess_list2, filter_list);

end