function [ts,ys] = DE2(f,t0,tN,y0,y1,h)
  ts = t0:h:tN;                    % defining time vector
  ys = zeros(length(ts));          % defining solution vector
  ys(1) = y0;                      % placing initial point in solution vector
  ys(2) = y0 + y1*h;               % placing second point in solution vector

  for i=2:length(ts)-1
      dy = (ys(i) - ys(i-1))/h;                         %finding dy
      ys(i+1) = 2*ys(i)-ys(i-1)+(h^2)*f(ts(i+1),ys(i),dy); %solving for y(n+1)
  end
end