function A = A1d(h, J)
% Central finite difference of d^2/dx^2 for 
% meshsize h and J points with homogeneous 
% Dirichlet BC.
e = ones(J, 1);
A = spdiags([-e/h^2 2*e/h^2 -e/h^2], [-1 0 1], J, J);
end