function g = pexa(x, t, T, nu)
% Manufactured solution for p
g = nu*sin(pi*x).*(exp(-pi^2*t) - exp(-pi^2*T)*(1+pi^2*t)/(1+pi^2*T));
end