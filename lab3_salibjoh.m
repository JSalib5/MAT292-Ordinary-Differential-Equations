%% ODE Lab: Creating your own ODE solver in MATLAB
% In this lab, you will write your own ODE solver for the Improved Euler method 
% (also known as the Heun method), and compare its results to those of |ode45|.
% 
% You will also learn how to write a function in a separate m-file and execute 
% it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each part using 
% cell mode to see the results. Compare the output with the PDF, which was generated 
% from this m-file.
% 
% There are six (6) exercises in this lab that are to be handed in on the due 
% date. Write your solutions in the template, including appropriate descriptions 
% in each step. Save the .m files and submit them online on Quercus.
%% Student Information
% Student Name: John Salib
% 
% Student Number: 1006729729
%% Creating new functions using m-files.
% Create a new function in a separate m-file:
% 
% Specifics: Create a text file with the file name f.m with the following lines 
% of code (text):
%%
% 
%  function y = f(a,b,c) 
%  y = a+b+c;
%
%% 
% Now MATLAB can call the new function f (which simply accepts 3 numbers and 
% adds them together). To see how this works, type the following in the matlab 
% command window: sum = f(1,2,3)
%% Exercise 1
% Objective: Write your own ODE solver (using the Heun/Improved Euler Method).
% 
% Details: This m-file should be a function which accepts as variables (t0,tN,y0,h), 
% where t0 and tN are the start and end points of the interval on which to solve 
% the ODE, y0 is the initial condition of the ODE, and h is the stepsize. You 
% may also want to pass the function into the ODE the way |ode45| does (check 
% lab 2).
% 
% Note: you will need to use a loop to do this exercise. You will also need 
% to recall the Heun/Improved Euler algorithm learned in lectures. 

%IEuler.m is the created improved euler method function
%% Exercise 2
% Objective: Compare Heun with |ode45|.
% 
% Specifics: For the following initial-value problems (from lab 2, exercises 
% 1, 4-6), approximate the solutions with your function from exercise 1 (Improved 
% Euler Method). Plot the graphs of your Improved Euler Approximation with the 
% |ode45| approximation.
% 
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
% 
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
% 
% (c) |y' = 1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
% 
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
% 
% Comment on any major differences, or the lack thereof. You do not need to 
% reproduce all the code here. Simply make note of any differences for each of 
% the four IVPs.

%a)
f = @(t,y) y*tan(t) + sin(t);
y0 = -0.5;
t0 = 0;
t1 = pi;
h = 0.001;
t = t0:h:t1;

soln = ode45(f, [t0, t1], y0);
my_soln = IEuler(f,t0,t1,y0,h);
plot(t, my_soln, soln.x, soln.y)

xlabel('t');
ylabel('y(t)');
legend('ode45', 'IEM','Location','Best');
snapnow
%graph of IEM method closely replicates that of the ODE45 function

%b)
f = @(t,y) 1/y^2;
y0 = 1;
t0 = 1;
t1 = 10;
h = 0.001;
t = t0:h:t1;

soln = ode45(f, [t0, t1], y0);
my_soln = IEuler(f,t0,t1,y0,h);
plot(t, my_soln, soln.x, soln.y)

xlabel('t');
ylabel('y(t)');
legend('ode45', 'IEM','Location','Best');
snapnow
%graph of IEM method closely replicates that of the ODE45 function

%c)
f = @(t,y) 1 - t*y/2;
y0 = -1;
t0 = 0;
t1 = 10;
h = 0.001;
t = t0:h:t1;

soln = ode45(f, [t0, t1], y0);
my_soln = IEuler(f,t0,t1,y0,h);
plot(t, my_soln, soln.x, soln.y)

xlabel('t');
ylabel('y(t)');
legend('ode45', 'IEM','Location','Best');
snapnow 
%graph of IEM method closely replicates that of the ODE45 function except
%for the reigon near t=1 to t=3 where the function reaches a max, at this
%point ODE and IEM depart from each other

%d)
f = @(t,y) y^3 - t^2;
y0 = 1;
t0 = 0;
t1 = 1;
h = 0.0001;
t = t0:h:t1;

soln = ode45(f, [t0, t1], y0);
my_soln = IEuler(f,t0,t1,y0,h);
plot(t, my_soln, soln.x, soln.y)

xlabel('t');
ylabel('y(t)');
legend('ode45', 'IEM','Location','Best');
%graph of IEM method closely replicates that of the ODE45 function up to
%t=0.5 where ODE 45 attempts to model the divergence of y to inf whereas
%IEM fails to model an point past t = 0.5
%% Exercise 3
% Objective: Use Euler's method and verify an estimate for the global error.
% 
% Details: 
% 
% (a) Use Euler's method (you can use euler.m from iode) to solve the IVP
% 
% |y' = 2 t sqrt( 1 - y^2 ) , y(0) = 0|
% 
% from |t=0| to |t=0.5|.
% 
% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.
% 
% (c) Read the attached derivation of an estimate of the global error for Euler's 
% method. Type out the resulting bound for En here in a comment. Define each variable.
% 
% (d) Compute the error estimate for |t=0.5| and compare with the actual error.
% 
% (e) Change the time step and compare the new error estimate with the actual 
% error. Comment on how it confirms the order of Euler's method.

%a)
f = @(t,y) 2*t*sqrt(1-y^2);
y0 = 0;
t0 = 0;
t1 = 0.5;
h = 0.01;
t = t0:h:t1;

soln = euler(f, y0, t);
%euler method's solution to the ODE using euler.m from IODE

%b)
y = sin(t.^2); %exact solution found by hand
%exact y at t = 0.5 
y_1 = sin(0.25);
%EM y at t = 0.5
ye_1 = soln(length(t));

%c)
%The bound for En <= (1+M)*dt/2 * (exp(M*dt*n)-1)
%where M = max(df/dt) in the range of t and dt = tN - t0 and n is the stepsize
%taken to go from t0 to tN (length of vector t)

%d)
M = 2; %maximum value of df/dt, f at t = 0.5
n = h; %step size between t0 and tN
dt = t1 - t0; %change in t
E_approx = (1+M)*dt/2 * (exp(M*dt*n)-1)
E_exact = y_1 - ye_1
%the approximate error upper bound is larger than the exact error, overestimating the
%deviation of the em function, this is to be expected as it is the upper
%bound of the possible error of the function

%e)
h = 0.1;
t = t0:h:t1;
soln = euler(f, y0, t);
ye_1 = soln(length(t));
E_approx = (1+M)*dt/2 * (exp(M*dt*h)-1)
E_exact = y_1 - ye_1

h = 0.001;
t = t0:h:t1;
soln = euler(f, y0, t);
ye_1 = soln(length(t));
E_approx = (1+M)*dt/2 * (exp(M*dt*h)-1)
E_exact = y_1 - ye_1

h = 0.0001;
t = t0:h:t1;
soln = euler(f, y0, t);
ye_1 = soln(length(t));
E_approx = (1+M)*dt/2 * (exp(M*dt*h)-1)
E_exact = y_1 - ye_1
%from the observed outputs it can be seen that both the approximate and
%exact errors remain at the same proportionality and change in magnitude
%with equal quantity, this confirms the first order relationship of Euler's
%method as the error change is linear
%% Adaptive Step Size
% As mentioned in lab 2, the step size in |ode45| is adapted to a specific error 
% tolerance.
% 
% The idea of adaptive step size is to change the step size |h| to a smaller 
% number whenever the derivative of the solution changes quickly. This is done 
% by evaluating f(t,y) and checking how it changes from one iteration to the next.
%% Exercise 4
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
% 
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as in 
% exercise 1, where |h| is an initial step size. You may also want to pass the 
% function into the ODE the way |ode45| does.
% 
% Create an implementation of Euler's method by modifying your solution to exercise 
% 1. Change it to include the following:
% 
% (a) On each timestep, make two estimates of the value of the solution at the 
% end of the timestep: |Y| from one Euler step of size |h| and |Z| from two successive 
% Euler steps of size |h/2|. The difference in these two values is an estimate 
% for the error.
% 
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be successful 
% and set the new solution value to be |Z+D|. This value has local error |O(h^3)|. 
% If |abs(D)>=tol|, reject this step and repeat it with a new step size, from 
% (c).
% 
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
% 
% Comment on what the formula for updating the step size is attempting to achieve.

%new function file created named AEuler.m
%updating the step size is attempting to reduce the value enough such that
%it maintains the method's output value within the tolerence error else it
%updates the stepsize as much as needed to be under the tolerence error.
%% Exercise 5
% Objective: Compare Euler to your Adaptive Euler method.
% 
% Details: Consider the IVP from exercise 3.
% 
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75| with 
% |h=0.025|.
% 
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.
% 
% (c) Plot both approximations together with the exact solution.

%a)
f = @(t,y) 2*t*sqrt(1-y^2);
y0 = 0;
t0 = 0;
t1 = 0.75;
h = 0.025;
t = t0:h:t1;

euler_soln = euler(f, y0, t); %euler's method solution to IVP

%b) 
[t_ea, y_ea] = AEuler(f,t0,t1,y0,h); %adaptive euler's method solution to IVP

%c)
plot(t, euler_soln, t_ea, y_ea, t_ea, sin(t_ea.^2))
xlabel('t');
ylabel('y(t)');
legend('EM', 'AEM', 'Exact', 'Location','Best');
%% Exercise 6
% Objective: Problems with Numerical Methods.
% 
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is closer 
% to the actual solution (done in 3.b)? Explain why.
% 
% (b) Plot the exact solution (from exercise 3.b), the Euler's approximation 
% (from exercise 3.a) and the adaptive Euler's approximation (from exercise 5) 
% from |t=0| to |t=1.5|.
% 
% (c) Notice how the exact solution and the approximations become very different. 
% Why is that? Write your answer as a comment.

%a) The adaptive Euler's Method solution is significantly closer (overlaps
%without zooming in) compared to the standard Euler's method. This is
%because the adaptive method continues to decrease the step size and make
%smaller incraments, extending the solution vectors to maintain an upper bound
%for the error of 1e-8. This keeps the function much closer to the exact
%solution.

%b)
f = @(t,y) 2*t*sqrt(1-y^2);
y0 = 0;
t0 = 0;
t1 = 1.5;
h = 0.025;
t = t0:h:t1;

euler_soln = euler(f, y0, t); %euler's method solution to IVP
[t_ea, y_ea] = AEuler(f,t0,t1,y0,h); %adaptive euler's method solution to IVP
plot(t, euler_soln, t_ea, y_ea, t_ea, sin(t_ea.^2))
xlabel('t');
ylabel('y(t)');
legend('EM', 'AEM', 'Exact', 'Location','Best');

%c) 
%past t = 1 the approximations and the exact solution diverge, this is
%beacuse the computation does not account for the imaginary components that
%appear past t = 1 in the derivative as it is restricted to 1 - y^2 > 0 for
%only real components in which both exact and approximations align closely.