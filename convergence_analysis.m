%Example template for analysis function
%INPUTS:
%solver_flag: an integer from 1-4 indicating which solver to use
% 1->Bisection 2-> Newton 3->Secant 4->fzero
%fun: the mathematical function that we are using the
% solver to compute the root of
%x_guess0: the initial guess used to compute x_root
%guess_list1: a list of initial guesses for each trial
%guess_list2: a second list of initial guesses for each trial
% if guess_list2 is not needed, then set to zero in input
%filter_list: a list of constants used to filter the collected data
function [root, p, k] = convergence_analysis(solver_flag, fun, ...
    x_guess0, guess_list1, guess_list2, filter_list)

    e_n0 = [];
    e_n1 = [];
    guess_it = [];

    % Compute the root once using a reliable method like fzero
    x_root = fzero(fun, x_guess0);
    fprintf('The reference root found is: %.14f\n', x_root);
    
    num_trials = length(guess_list1);
    for i = 1:num_trials
        x0_trial = [];
        x1_trial = [];
        it_trial = [];

        % Check if the solver_flag is valid and call the appropriate solver
        switch solver_flag
            case 1 % Bisection
                [~, x0_trial, x1_trial, it_trial] = bisection_solver(fun, guess_list1(i), guess_list2(i));
            case 2 % Newton
                % Create separate function handles for f and df/dx for the newton_solver
                f_handle = @(x) get_f(fun, x);
                dfdx_handle = @(x) get_dfdx(fun, x);
                [~, guesses] = newton_solver(f_handle, dfdx_handle, guess_list1(i), 1e-14, 1e-14, 100);
                if length(guesses) > 1
                    x0_trial = guesses(1:end-1);
                    x1_trial = guesses(2:end);
                    it_trial = 1:length(x0_trial);
                end
            case 3 % Secant
                [x1_trial, x0_trial, it_trial, ~] = secant_method(fun, guess_list1(i), guess_list2(i), 1e-14, 1e-14, 100);
            case 4 % fzero
                % NEEDS WORK
                root = fzero(fun, guess_list1(i));
                disp('NEEDS WORK.');
            otherwise
                error('Invalid solver_flag. Must be an integer from 1 to 4.');
        end
    
        % Append data from the current trial
        if ~isempty(x0_trial)
            e_n0 = [e_n0 abs(x0_trial - x_root)];
            e_n1 = [e_n1 abs(x1_trial - x_root)];
            guess_it = [guess_it it_trial];
        end
    end

   
    %example for how to filter the error data
    %currently have error_list0, error_list1, index_list
    %data points to be used in the regression
    x_regression = []; % e_n
    y_regression = []; % e_{n+1}
    %iterate through the collected data
    for n=1:length(guess_it)
         %if the error is not too big or too small
         %and it was enough iterations into the trial...
         if e_n0(n)>filter_list(1) && e_n0(n)<filter_list(2) && ...
            e_n1(n)>filter_list(3) && e_n1(n)<filter_list(4) && ...
            guess_it(n)>filter_list(5)
            %then add it to the set of points for regression
            x_regression(end+1) = e_n0(n);
            y_regression(end+1) = e_n1(n);
        end
    end

    [p, k] = generate_error_fit(x_regression, y_regression);

    loglog(e_n0,e_n1,'ro','markerfacecolor','r','markersize',1); %loglog sanity check
    hold on;
    %example for how to plot fit line
    %generate x data on a logarithmic range
    fit_line_x = 10.^(-16:.01:1);
    %compute the corresponding y values
    fit_line_y = k*fit_line_x.^p;
    %plot on a loglog plot.
    loglog(x_regression, y_regression,'bo','markerfacecolor','b','markersize',1);
    loglog(fit_line_x,fit_line_y,'k-','linewidth',2)

    title('Convergence Analysis');
    xlabel('\epsilon_n');
    ylabel('\epsilon_{n+1}');
    legend('Raw Data', 'Regression Data', 'Filtered Data', 'Location', 'NorthWest');
    hold off;
    fprintf('Regression coefficients are:\n');
    fprintf('p: %.3f\n', p);
    fprintf('k: %.3f\n', k);

    root = x_root;
end

% Helper function to get only the function value
function f = get_f(fun_handle, x)
    [f, ~] = fun_handle(x);
end

% Helper function to get only the derivative value
function dfdx = get_dfdx(fun_handle, x)
    [~, dfdx] = fun_handle(x);
end

function [f, dfdx] = fun(x)
    f =  (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) - 0.7 - exp(x/6);
    dfdx =  (3*x.^2)/100 - (x)/4 + 2 + 3*cos(x/2 + 6)*0.5 - (1/6)*exp(x/6);
end
%example for how to compute the fit line
%data points to be used in the regression
%x_regression -> e_n
%y_regression -> e_{n+1}
%p and k are the output coefficients
function [p,k] = generate_error_fit(x_regression,y_regression)
    %generate Y, X1, and X2
    %note that I use the transpose operator (')
    %to convert the result from a row vector to a column
    %If you are copy-pasting, the ' character may not work correctly
    Y = log(y_regression)';
    X1 = log(x_regression)';
    X2 = ones(length(X1),1);
    %run the regression
    coeff_vec = regress(Y,[X1,X2]);
    %pull out the coefficients from the fit
    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end
%example of how to implement finite difference approximation
%for the first and second derivative of a function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%x: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%dfdx: approximation of fun'(x)
%d2fdx2: approximation of fun''(x)
function [dfdx,d2fdx2] = approximate_derivative(fun,x)
    %set the step size to be tiny
    delta_x = 1e-6;
    %compute the function at different points near x
    f_left = fun(x-delta_x);
    f_0 = fun(x);
    f_right = fun(x+delta_x);
    %approximate the first derivative
    dfdx = (f_right-f_left)/(2*delta_x);
    %approximate the second derivative
    d2fdx2 = (f_right-2*f_0+f_left)/(delta_x^2);
end