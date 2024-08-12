function [fixedPoint, err, iterations] = fixed_point_iteration(g, x0, TOL, max_iter)
   % g: The function handle for which we want to find the fixed point.
   % x0: Initial guess for the fixed point.
   % TOL: Tolerance (stop when |x_n+1 - x_n| < tol).
   % max_iter: Maximum number of iterations.
   % Check for the number of input arguments
   n = nargin;
   if n < 2
       error("This function must have at least 2 inputs");
   end
   % Set default values if not provided
   if n < 3 || isempty(TOL)
       TOL = 2e-10;
   end
   if n < 4 || isempty(max_iter)
       max_iter = 10000;
   end
   %Initialize variables
   x = x0;
   iter = 0;
   err = inf; 
   while err > TOL && iter < max_iter
       x_i = g(x); 
       err = abs(g(x_i) - g(x));
       x = x_i;
       iter = iter + 1;
   end
   % If max_iter iterations are reached without convergence, return the last estimate.
   fixedPoint = x;
   iterations = iter;
end
