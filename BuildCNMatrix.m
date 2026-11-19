function M = BuildCNMatrix(A, J, dt, N, nu, Ay0, Ap0, AyT, ApT)
% Build the matrix M to solve the problem MU = F,
% with U = [Y; P] for both heat and 
% advection-diffusion constraints.
% The time integration scheme is Crank-Nicolson. 
% There are four components in the matrix 
% M = [M11, M12; M21, M22]. The first block of 
% M11 and M12 and the last block of M21 and M22
% are specially adapted to different transmission
% conditions. 
%
% Input
% ----------
% A : matrix
%     (J-1)*(J-1) space matrix
% J : int
%     number of points in space (J+1)
% dt : float
%     time step
% N : int
%     number of points in time (N+1)
% nu : float
%     penalization for control
% Ay0 : matrix
%     (J-1)*(J-1) matrix -- first block of M11
% Ap0 : matrix
%     (J-1)*(J-1) matrix -- first block of M12
% AyT : matrix
%     (J-1)*(J-1) matrix -- last block of M21
% ApT : matrix
%     (J-1)*(J-1) matrix -- last block of M22

e = ones(N+1, 1);

% construct four matrices
M11 = kron(speye(N+1), speye(J-1)+dt*A/2) - kron(spdiags(e, -1, N+1, N+1), speye(J-1)-dt*A/2);
M22 = -M11';
M21 = kron(spdiags([e e], [0 1], N+1, N+1), -dt*speye(J-1)/2);
M12 = M21'/nu;

% adapt the initial and final condition
M11(1:J-1, 1:J-1) = Ay0;
M12(1:J-1, 1:J-1) = Ap0;
M21(end-J+2:end, end-J+2:end) = AyT;
M22(end-J+2:end, end-J+2:end) = ApT;

% concatenation of four matrices
M = [M11, M12; M21, M22];
end