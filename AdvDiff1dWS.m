% --------------------------------------------- %
% Script to generate data of the weak scalability
% test for a 1D advection-diffusion OCP.
% --------------------------------------------- %
clc; clear; close all

nu = 1e-1;
c = 1;
kap = 1e-2;

L = 1; 
hx = 2^-5;
x = 0:hx:L; 
Nx = length(x) - 1;

ht = hx;

Pe = c*L/kap %Peclet number

Twl = 2^-1; 
Ntw = 2.^(1:1:9); 

Nit = 20;  
erry = ones(length(Ntw), Nit); errp = erry; 

filename = sprintf('AdvDiff1dWSdata2-512_nu_1e%d_c_%d_kap_1e%d_L_%d_hxht_2^%d_Twl_2^%d_Nit_%d.mat',log10(nu), c, log10(kap), L, log2(hx), log2(Twl), Nit);

if isfile(filename)
  disp('Data exists and script stopped.');
  return;
end

A = A1dAdv(c, kap, hx, Nx-1);
rng(0);

for k = 1:length(Ntw)
  Nt = Ntw(k)*Twl/ht

  Mm = BuildCNMatrix(A, Nx, ht, Twl/ht, nu, speye(Nx-1), sparse(Nx-1, Nx-1), sparse(Nx-1, Nx-1), speye(Nx-1)); 
  dMm = decomposition(Mm, 'lu');

  gl = rand(Nx-1, Ntw(k)-1); gr = gl;
  for nit = 1:Nit
    nit
    gltp = zeros(size(gl)); grtp = gltp;
    ytp = []; ptp = [];
    for m = 1:Ntw(k)
      if m == 1
        F1 = BuildCNRHS(Nx, ht, Twl/ht, sparse(Twl/ht+1,Nx-1), sparse(Nx-1, 1), gr(:, 1));
        X1 = dMm\F1; X1 = reshape(X1, Nx-1, 2*Twl/ht+2);
        y1 = X1(:, 1:Twl/ht+1); p1 = X1(:, Twl/ht+2:end);

        ytp = [ytp y1]; ptp = [ptp, p1(:, 1:end-1)];
        gltp(:, 1) = y1(:, end);

      elseif m == Ntw(k)
        Fn = BuildCNRHS(Nx, ht, Twl/ht, sparse(Twl/ht+1,Nx-1), gl(:, end), sparse(Nx-1, 1));
        Xn = dMm\Fn; Xn = reshape(Xn, Nx-1, 2*Twl/ht+2);
        yn = Xn(:, 1:Twl/ht+1); pn = Xn(:, Twl/ht+2:end);

        ytp = [ytp yn(:,2:end)]; ptp = [ptp, pn];
        grtp(:, end) = pn(:, 1);

      else
        Fm = BuildCNRHS(Nx, ht, Twl/ht, sparse(Twl/ht+1,Nx-1), gl(:, m-1), gr(:, m));
        Xm = dMm\Fm; Xm = reshape(Xm, Nx-1, 2*Twl/ht+2);
        ym = Xm(:, 1:Twl/ht+1); pm = Xm(:, Twl/ht+2:end);

        ytp = [ytp ym(:,2:end)]; ptp = [ptp, pm(:, 1:end-1)];
        gltp(:, m) = ym(:, end); grtp(:, m-1) = pm(:, 1);
      end
    end

    gl = gltp; gr = grtp;
    erry(k,nit) = norm(ht*full(ytp));
    errp(k,nit) = norm(ht*full(ptp));
  end
end

save(filename, 'nu', 'c', 'kap', 'L', 'hx', 'x', 'Nx', 'ht', 'Twl', 'Ntw', 'Nit', 'erry', 'errp', 'A');