function main()

    %first root
    x_left = -15;
    x_right = 10;
    
    [bi_root, bi_c, bi_fc, bi_guess_it] = bisection_solver(@fun, x_left, x_right);
    bi_n_current; %current c guesses for bisection
    bi_n_next; %next c guesses for bisection
    bi_guess_it; %iteration number for bisection
    
    fprintf('The first root found is: %.14f\n', bi_root);
    
    % %second root
    % x_left = 25;
    % x_right = 40;
    % root2 = bisection_solver(@fun, x_left, x_right);
    % fprintf('The second root found is: %.14f\n', root2);

    x_root = bi_root; % setting true root as root found via bisection, can replace with any other method if desired. fzero yields exact same precision result

end

function [f, dfdx] = fun(x)
    f =  (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) - 0.7 - exp(x/6);
    dfdx =  (3*x.^2)/100 - (x)/4 + 2 + 3*cos(x/2 + 6)*0.5 - (1/6)*exp(x/6);
end
