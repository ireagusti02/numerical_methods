function [sol, err, iter] = bisection(f, a, b, TOL, maxIter)
%f should be a function handle which we want the numerical solution of.
% a and b are the bounds of the interval in which we want the solution in
%tol is the tolerated error we are willing to take
%maxIter is the max amout of times we are willing to run this.
%sol is the solution
%err is the error after iter, which is the number of iterations.
   % Check for the number of input arguments
   n = nargin;
   if n < 3
       error("This function must have at least 3 inputs");
   end
   % Set default values if not provided
   if n < 4
       TOL = 2e-10;
   end
   if n < 5
       maxIter = 100;
   end
   % Initialize variables
   err = inf;
   sol = NaN;
   iter = 0;
   % Check if the interval [a, b] is valid
   if f(a) * f(b) > 0
       fprintf('There is no guaranteed root in this interval\n');
       return;
   end
   % Bisection method
   while err > TOL && iter < maxIter
       c = (a + b) / 2;
       product = f(a) * f(c);
       if product < 0
           b = c;
       elseif product > 0
           a = c;
       else
           sol = c;
           err = 0;
       end
       iter = iter + 1;
       err = abs(f(sol));
   end
   % Set the solution if found within the tolerance
   if err <= TOL
       sol = (a + b) / 2;
   end
