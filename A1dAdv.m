function A = A1dAdv(c, kap, h, J)
% Finite difference approximation of 
% c*d/dx - kap*d^2/dx^2 with c>0 and kap>0 for 
% meshsize h and J points with homogeneous 
% Dirichlet BC. Central finite difference for
% d^2/dx^2 and upwind for d/dx.
e = ones(J, 1);
A = spdiags([(-kap/h^2 - c/h)*e, (2*kap/h^2 + c/h)*e, -kap/h^2*e], [-1 0 1], J, J);
end
