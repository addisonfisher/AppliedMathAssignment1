%Example sigmoid function
function [f_val,dfdx] = test_function03(x)
    global input_list;
    input_list(:,end+1) = x;
    a = 27.3; b = 2; c = 8.3; d = -3;
    H = exp((x-a)/b);
    dH = H/b;
    L = 1+H;
    dL = dH;
    f_val = c*H./L+d;
    dfdx = c*(L.*dH-H.*dL)./(L^2);

    x = 0:0.1:50;
    y = @(x) c*exp((x-a)/b)./(1+exp((x-a)/b)) + d;
    
    num_trials = 200;
    guess_list1 = randi([-5, 5], 1, num_trials);
    guess_list2 = randi([45, 50], 1, num_trials);

    filter_list = [1e-15, 1e-2, 1e-14, 1e-2, 2];
    root = convergence_analysis(1, y, 5, guess_list1, guess_list2, filter_list);

    figure();
    plot(x, y(x), 'LineWidth', 2);
    hold on;
    plot(root, y(root), 'MarkerSize', 4, 'Marker', 'o', 'MarkerFaceColor', 'r');
    y_value = y(root);
    yline(y_value, 'k--', 'LineWidth', 1.5);

end