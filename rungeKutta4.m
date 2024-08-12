function [yout, approxList] = rungeKutta4(f, t0, y0, tn, h)
    % f: Function handle for the system of ODEs (dy/dt = f(t, y))
    % t0: Initial time
    % y0: Initial condition
    % tn: End time
    % h: Step size

    % Check for correct input
    [r, c] = size(y0);
    if ((r > 1) && (c > 1))
        error("input is a matrix but should be a vector")
    elseif (r > 1)
        y0 = y0';
    end

    % Check for correct output
    [r1, c1] = size(f(t0, y0));
    if ((r1 > 1) && (c1 > 1))
        error("output is a matrix but should be a vector")
    elseif (r1 > 1)
        f = @(t, y) f(t, y)';
    end

    % Initialize variables
    n = ceil((tn - t0) / h);
    t_values = linspace(t0, tn, n + 1);
    steps = length(t_values);
    step = t_values(2) - t_values(1);
    half_step = step / 2;
    y = y0;

    % Initialize approxList
    approxList = zeros(steps, length(y0) + 1);
    approxList(1, 1) = t0;
    approxList(1, 2:end) = y0;

    % Runge-Kutta 4th order integration (vectorized)
    for i = 2:steps
        t = t_values(i);

        % Calculate the four stages
        k1 = step * f(t, y);
        k2 = step * f(t + half_step, y + k1 / 2);
        k3 = step * f(t + half_step, y + k2 / 2);
        k4 = step * f(t + step, y + k3);

        % Update y using the weighted sum of the stages
        y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;

        % Store the current iteration in approxList
        approxList(i, :) = [t, y];
     
    end

    % Output the approximation at tn
    yout = approxList(end, 2:end);
end
