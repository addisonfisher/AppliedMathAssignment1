function egg()
    % %set the oval hyper-parameters
    egg_params = struct();
    egg_params.a = 3;
    egg_params.b = 2;
    egg_params.c = .15;
    % 
    %specify the position and orientation of the egg
    x0 = 5; y0 = 5; theta = pi/6;

    %set up the axis
    hold on;
    axis equal;
    axis square;
    axis([0, 10, 1, 10]);
    s_perimeter = linspace(0, 1, 100);
    % 
    % 
    [V_vals, G_vals] = egg_func(s_perimeter, x0, y0, theta, egg_params);
    plot(V_vals(1,:), V_vals(2,:), '-', 'LineWidth', 1);
    %plot(V_vals(1,:), V_vals(2,:), 'ro', 'markerfacecolor', 'r', 'markersize', 4);
    % 
    % s_tangent = .3;
    % [V_tangent, G_tangent] = egg_func(s_tangent, x0, y0, theta, egg_params)
    % 
    % plot(V_tangent(1), V_tangent(2), 'ro', 'markerfacecolor', 'r', 'markersize', 4);
    % plot(V_tangent(1) + [0, G_tangent(1)], V_tangent(2) + [0, G_tangent(2)], 'k',  'LineWidth', 2);    
    % 
    % [xmin, xmax, ymin, ymax] = find_bounding_box(x0, y0, theta, egg_params)
    % 
    % plot ([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], 'b-', 'LineWidth',2);
    % hold off;
    % 
    % y_ground = 1;
    % x_wall = 15;
   
    animation_example(@egg_trajectory01, 0, 25);

end

function[xmin, xmax, ymin, ymax] = find_bounding_box(x0, y0, theta, egg_params)
    egg_wrapper_func2_x = @(s_in) egg_wrapper_func1_x(s_in, x0, y0, theta, egg_params);
    egg_wrapper_func2_y = @(s_in) egg_wrapper_func1_y(s_in, x0, y0, theta, egg_params);

    s_guess_list = 0:0.1:1;

    dxtol = 1e-14; ytol = 1e-14; max_iter = 200; dfdxmin = 1e-8;

    x_list = [];
    y_list = [];
    for s_guess = s_guess_list
            s_rootx = secant_method(egg_wrapper_func2_x, s_guess, s_guess + 1e-4);%, dxtol, ytol, max_iter, dfdxmin);
            [V, G] = egg_func(s_rootx, x0, y0, theta, egg_params);
            x_list(end + 1) = V(1);

            s_rooty = secant_method(egg_wrapper_func2_y, s_guess, s_guess+1e-4); %, dxtol, ytol, max_iter, dfdxmin);
            [V, G] = egg_func(s_rooty, x0, y0, theta, egg_params);
            y_list(end + 1) = V(2);
    end
    xmin = min(x_list);
    xmax = max(x_list);
    ymin = min(y_list);
    ymax = max(y_list);

end

%This first wrapping is going to make sure the output is the single scalar,
%x component of G
function Gx = egg_wrapper_func1_x(s, x0, y0, theta, egg_params)
    [~, G] = egg_func(s, x0, y0, theta, egg_params);
    Gx = G(1);
end

function Gy = egg_wrapper_func1_y(s, x0, y0, theta, egg_params)
    [~, G] = egg_func(s, x0, y0, theta, egg_params);
    Gy = G(2);
end

%This function generates the parametric curve describing an oval
%INPUTS:
%s: the curve parametr. s is a number from 0 to 1. The curve function has a
% period of 1, so s=.3 and s=1.3 will generate the same output
% s can also be a list (row vector) of numbers
%theta: rotation of the oval. theta is a number from 0 to 2*pi.
% Increasing theta rotates the oval counterclockwise
%x0: horizontal offset of the oval
%y0: vertical offset of the oval
%egg_params: a struct describing the hyperparameters of the oval
% egg_params has three variables, a,b, and c
% without any rotation/translation, the oval satisfies the equation:
% xˆ2/aˆ2 + (yˆ2/bˆ2)*eˆ(c*x) = 1
% tweaking a,b,c changes the shape of the oval
%OUTPUTS:
%V: the position of the point on the oval given the inputs
% If s is a single number, then V will have the form of a column vector
% [x_out;y_out] where (x_out,y_out) are the coordinates of the point on
% the oval. If the input t is a list of numbers (a row vector) i.e.:
% s = [s_1,...,s_N]
% then V will be an 2xN matrix:
% [x_1,...,x_N; y_1,...,y_N]
% where (x_i,y_i) correspond to input s_i
%G: the gradient of V taken with respect to s
% If s is a single number, then G will be the column vector [dx/ds; dy/ds]
% If s is the list [s_1,...,s_N], then G will be the 2xN matrix:
% [dx_1/ds_1,...,dx_N/ds_N; dy_1/ds_1,...,dy_N/ds_N]
function [V, G] = egg_func(s,x0,y0,theta,egg_params)
    %unpack the struct
    a=egg_params.a;
    b=egg_params.b;
    c=egg_params.c;
    %compute x (without rotation or translation)
    x = a*cos(2*pi*s);
    %useful intermediate variable
    f = exp(-c*x/2);
    %compute y (without rotation or translation)
    y = b*sin(2*pi*s).*f;
    %compute the derivatives of x and y (without rotation or translation)
    dx = -2*pi*a*sin(2*pi*s);
    df = (-c/2)*f.*dx;
    dy = 2*pi*b*cos(2*pi*s).*f + b*sin(2*pi*s).*df;
    %rotation matrix corresponding to theta
    R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    %compute position and gradient for rotated + translated oval
    V = R*[x;y]+[x0*ones(1,length(theta));y0*ones(1,length(theta))];
    G = R*[dx;dy];
end

function [x0,y0,theta] = egg_trajectory01(t)
    x0 = 7*t + 8;
    y0 = -6*t.^2 + 20*t + 6;
    theta = 5*t;
end

function ymin = find_ymin(traj_func, t, egg_params)
    [x0, y0, theta] = traj_func(t);

    s_perimeter = linspace(0, 1, 1000);
    [V, ~] = egg_func(s_perimeter, x0, y0, theta, egg_params);

    ymin = min(V(2, :));

end

function xmax = find_xmax(traj_func, t, egg_params)
    [x0, y0, theta] = traj_func(t);

    s_perimeter = linspace(0, 1, 1000);
    [V, ~] = egg_func(s_perimeter, x0, y0, theta, egg_params);

    xmax = max(V(1, :));
end


%Function that computes the collision time for a thrown egg
%INPUTS:
%traj_fun: a function that describes the [x,y,theta] trajectory
% of the egg (takes time t as input)
%egg_params: a struct describing the hyperparameters of the oval
%y_ground: height of the ground
%x_wall: position of the wall
%OUTPUTS:
%t_ground: time that the egg would hit the ground
%t_wall: time that the egg would hit the wall
function [t_ground, t_wall] = collision_func(traj_fun, egg_params, y_ground, x_wall)
    %Pass in the function of x_wall - x_max (which we take from the
    %bounding box)
        
    fy_min = @(t) find_ymin(traj_fun, t, egg_params) - y_ground;
    fx_max = @(t) find_xmax(traj_fun, t, egg_params) - x_wall;
    
    % initial guess
    x_guess0 = 10;
    
    % guess lists
    num_trials = 250;
    tspan = linspace(0, 20, 250);   % search from t=0 to 20

    guess_list1 = tspan(1:end-1);   % left guesses
    guess_list2 = tspan(2:end);     % right guesses
    
    % filtering constants
    filter_list = [1e-15, 1e-2, 1e-14, 1e-2, 2];
    solver_flag = 2;
    t_ground = convergence_analysis(solver_flag, fy_min, x_guess0, guess_list1, guess_list2, filter_list);
    t_wall = convergence_analysis(solver_flag, fx_max, x_guess0, guess_list1, guess_list2, filter_list);

end

%Short example demonstrating how to create a MATLAB animation
%In this case, a square moving along an elliptical path
function animation_example(traj_func, y_ground, x_wall)
%Define the coordinates of the square vertices (in its own frame)

%set up the plotting axis
hold on; axis equal; axis square
axis([0,40,0,40])
%initialize the plot of the square
egg_plot = plot(0,0,'r');
box_plot = plot(0,0,'k');

wall_plot = plot(0,0,'k', LineWidth=2);
ground_plot = plot(0,0,'k', LineWidth=2);

%Initialize egg params
egg_params = struct();
egg_params.a = 3;
egg_params.b = 2;
egg_params.c = .15;

[t_ground, t_wall] = collision_func(traj_func, egg_params, y_ground, x_wall);

if t_ground < 0 
    t_ground = inf;
end

if t_wall < 0 
    t_wall = inf;
end

t_stop = min(t_ground, t_wall); 
hold on;

%iterate through time
for t=0:.001:t_stop %abs(time)
    [position_x, position_y, theta] = traj_func(t);
    s_perimeter = linspace(0, 1, 100);
    [V_vals, ~] = egg_func(s_perimeter, position_x, position_y, theta, egg_params);

    %compute positions of square vertices (in world frame)
    x_plot = V_vals(1,:);
    y_plot = V_vals(2,:);
    
    %update the coordinates of the square plot
    set(egg_plot,'xdata',x_plot,'ydata',y_plot);
    
    [xmin, xmax, ymin, ymax] = find_bounding_box(position_x, position_y, theta, egg_params);

    
    set(box_plot,'xdata',[xmin, xmax, xmax, xmin, xmin],'ydata', [ymin, ymin, ymax, ymax, ymin]);
    %update the actual plotting window

    
    set(wall_plot,'xdata',[x_wall, x_wall],'ydata', [0, 40]);

    set(ground_plot,'xdata',[0, 40],'ydata', [y_ground, y_ground]);

    drawnow;
end

end