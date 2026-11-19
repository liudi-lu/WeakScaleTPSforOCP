function g = yexa(x, t, T)
% Manufactured solution for y
g = sin(pi*x).*(t.*exp(-pi^2*t) - exp(-pi^2*T)*t/(1+pi^2*T));
end