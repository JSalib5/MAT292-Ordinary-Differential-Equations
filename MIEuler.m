function [ts,ys] = IEuler(f,g,t0,tN,x0,h)
  ts = t0:h:tN;                    % defining time vector
  ys = zeros(2,length(ts));        % defining solution matrix
  ys(1,1) = x0(1);                 % placing initial point in solution matrix
  ys(2,1) = x0(2);

  for i=1:length(ts)-1
    s1 = f(ts(i),ys(1,i),ys(2,i));         % evaluate ODE1 at current point
    h1 = g(ts(i),ys(1,i),ys(2,i));         % evaluate ODE1 at current point
    yE1 = ys(1,i) + s1*h;                  % find Euler Method value for ODE 1
    yE2 = ys(2,i) + h1*h;                  % find Euler Method value for ODE 2
    s2 = f(ts(i)+h,yE1,yE2);               % evaluate ODE 1 at Euler Method result
    h2 = g(ts(i)+h,yE1,yE2);               % evaluate ODE 2 at Euler Method result
    ys(1,i+1) = ys(1,i) + (s1+s2)/2*h;     % use mean of both points (s1+s2)/2 to find IEM at t and storing in solution
    ys(2,i+1) = ys(2,i) + (h1+h2)/2*h;
  end
end