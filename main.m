% --- Main Script ---

test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) - 0.7 - exp(x/6);

%first root
left_bound = -15;
right_bound = 10;

a_thresh = 1e-14;
b_thresh = 1e-14;
max_it = 1000;

root1 = bisection_method(test_func01, left_bound, right_bound, a_thresh, b_thresh, max_it);

fprintf('The first root found is: %.14f\n', root1);

%second root
left_bound = 25;
right_bound = 40;
root2 = bisection_method(test_func01, left_bound, right_bound, a_thresh, b_thresh, max_it);
fprintf('The second root found is: %.14f\n', root2);