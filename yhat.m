function g = yhat(x, t, T, nu)
% Target function for manufactured solutions
g = nu*sin(pi*x).*((t/nu+2*pi^2).*exp(-pi^2*t) - (t/nu+pi^4*t)*exp(-pi^2*T)/(1+pi^2*T));
end