function [rhoT, normT, boundT] = Computebound(d, nu, Twl, Ntw)
%define quantities in the paper
sgm = sqrt(d^2 + 1/nu); 
sh = @(tj) sinh(sgm * tj);  
ch = @(tj) cosh(sgm * tj);  
C1dt = @(dt)-(1/nu)*sh(dt)/(sgm*ch(dt)+d*sh(dt)); 
C2dt = @(dt)sgm/(sgm*ch(dt)+d*sh(dt)); 

% 2-by-2 block matrices T_{l,m}, T_{d,m}, T_{r,m}
Tl = [0 0; 0 C2dt(Twl)];
Td = [0 C1dt(Twl); -nu*C1dt(Twl) 0];
Tr = [C2dt(Twl) 0; 0 0];

rhoT = []; normT = [];
for N = 2:Ntw
  e = ones(N, 1);

  %define matrix T^{PS}_{N,m}
  T = kron(spdiags(e, -1, N, N), Tl) + kron(speye(N), Td) + kron(spdiags(e, 1, N, N), Tr);

  %spectral radius
  rhoT = [rhoT max(abs(eig(full(T))))];
  
  %infinity norm
  normT = [normT norm(T, 'inf')];      
end

%our bound
boundT = sqrt(nu*C1dt(Twl)^2+C2dt(Twl)^2)*ones(size(rhoT));
end