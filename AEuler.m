function [ts,ys] = AEuler(f,t0,tN,y0,h)
  ys(1) = y0;                       % placing initial point in solution vector
  ts(1) = t0;
  tol = 1e-8;                      % assigning max tolerence error
  
  while ts(end) < tN               % iterating through until end of time interval is reached
    s1 = f(ts(end),ys(end));       % evaluate ODE at current point
    Y = ys(end) + s1*h;            % finding the Euler value at end of 1 h time step
    Z = ys(end) + s1*(h/2);
    s2 = f(ts(end)+h/2,Z);
    Z = Z + s2*(h/2);              % finidng Euler value after 2 successive steps of h/2
    D = Z-Y;                       % assigning error for current step
    
    if abs(D)<tol
        ys(end+1) = Z-D;           % successful step hence assigning new y value
        ts(end+1) = ts(end) + h;   % assigning new time value
    else 
        h = 0.9*h*min(max(tol/abs(D),0.3),2); % assigning new time step
        continue                   %skips storage of new result and restarts loop
    end
    
  end
end