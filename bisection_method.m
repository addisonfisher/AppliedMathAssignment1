% left_start = -15;
% right_start = 20;
% middle_start = (left_start + right_start)/2;
% 
% a = left_start;
% b = right_start;
% c = middle_start;
% 
% interval = check_signs(a, b, c);
% 
% function interval = check_signs(a, b, c)
% 
%     %Definition of the test function and its derivative
%     test_func01 = @(x) (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
%     test_derivative01 = @(x) 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
% 
%     fa = test_func01(a);
%     fb = test_func01(b);
%     fc = test_func01(c);
% 
%     if(sign(fa) == sign(fc))
%         interval = [c, b];
%     elseif(sign(fb) == sign(fc))
%         interval = [a, c];
%     else    
%         disp('error');
%         interval = [];
%     end
% end


function root = bisection_method(fun, a, b, a_thresh, b_thresh, max_it)
    if sign(fun(a)) == sign(fun(b))
        error('The initial interval does not cross zero. Bisection method will not work!');
    end
   
    c_old = a;

    for i = 1:max_it
        c = (a + b) / 2;
        
        delta_x = abs(c - c_old);
        f_c = fun(c);

        if delta_x < a_thresh || abs(f_c) < b_thresh
            root = c;
            fprintf('Converged after %d iterations.\n', i);
            return;
        end
        
        if sign(fun(a)) == sign(f_c) 
            b = c;
        else
            a = c;
        end
        
        c_old = c;
    end
end