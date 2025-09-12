function test

    % 1. Set the solver_flag for the bisection method.
    solver_flag = 3;
    
    % 2. Define the handle for the mathematical function to be analyzed.
    % This function is provided in the assignment PDF  and your convergence_analysis.m file.
    fun = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) - 0.7 - exp(x/6);
    
    % 3. Provide an initial guess to find the "true" root for error calculation.
    % From the plot in the PDF, we can see a root is near -5[cite: 32].
    x_guess0 = -5;
    
    % 4. Generate lists of initial guesses for the trials.
    % For the bisection method, these are the left and right interval endpoints.
    % We'll create 200 random intervals that should contain the root.
    num_trials = 200;
    guess_list1 = randi([-30, -5], 1, num_trials); % x_left values
    guess_list2 = randi([5, 30], 1, num_trials);    % x_right values
    
    % 5. Define the constants for filtering the error data, as shown in the PDF[cite: 248].
    % [e_n0_min, e_n0_max, e_n1_min, e_n1_max, min_iterations]
    filter_list = [1e-15, 1e-2, 1e-14, 1e-2, 2];
    
    % 6. Call the analysis function with the defined inputs.
    convergence_analysis(solver_flag, fun, x_guess0, guess_list1, guess_list2, filter_list);

end