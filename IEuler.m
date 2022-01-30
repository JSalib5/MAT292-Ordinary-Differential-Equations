function ys = IEuler(f,t0,tN,y0,h)
  ys = t0:h:tN;                    %defining solution vector
  t = t0; 
  y = y0;                          % defining initial point
  ys(1) = y;                       %placing initial point in solution vector

  for i=1:length(ys)-1
    s1 = f(t,y);                   % evaluate ODE at current point
    yE = y + s1*h;                 % find Euler Method value
    s2 = f(t+h,yE);                % evaluate ODE at Euler Method result
    y = y + (s1+s2)/2*h;           % use mean of both points (s1+s2)/2 to find IEM at t
    t = t + h;                     % cycling t
    ys(i+1) = y;                   % store result in vector ys
  end
end