function y = RK4_new(F, tspan, y0, h, params)
% F: Function handle for the derivative function
% tspan: Time span [t_start, t_end]
% y0: Initial value
% h: Step size

t_start = tspan(1);
t_end = tspan(end);
n = round((t_end - t_start) / h); % Number of steps
t = t_start; % Initial time
y = zeros(length(y0), n+1); % Array to store solution
y(:,1) = y0; % Initial condition

for i = 1:n
    k1 = F(t, y(:,i), params);
    k2 = F(t + h/2, y(:,i) + h/2 * k1, params);
    k3 = F(t + h/2, y(:,i) + h/2 * k2, params);
    k4 = F(t + h, y(:,i) + h * k3, params);
    y(:,i+1) = y(:,i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    t = t + h;
end