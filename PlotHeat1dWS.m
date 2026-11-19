% --------------------------------------------- %
% Script to generate the figure of the weak
% scalability test for a 1D heat OCP with 
% manufactured solutions, and the comparaison 
% of spectral radius, infinity norm with our   
% bound. The time window length (Twl) needs to be
% change to generate different figures. When
% saving the figure, change also Twl value in
% the filename manually!
% --------------------------------------------- %
clc; clear; close all

nu = 1e-1;

L = 1; 
hx = 2^-5;
x = 0:hx:L;
Nx = length(x) - 1;
      
ht = hx;

Twl = 2^-1; %2^-2;2^-3;2^-4;
Ntw = 2.^(1:1:4);

A = A1d(hx, Nx-1);
dvec = eig(A);                
d = min(abs(dvec));
[rhoT, normT, boundT] = Computebound(d, nu, Twl, Ntw(end));

Nit = 10;                   
erry = ones(length(Ntw), Nit); errp = erry;

for k = 1:length(Ntw)
  %n = Ntw(k);
  t = 0:ht:Ntw(k)*Twl;
  Nt = length(t) - 1;

  %solve for one domain reference solutions
  M = BuildCNMatrix(A, Nx, ht, Nt, nu, speye(Nx-1), sparse(Nx-1, Nx-1), sparse(Nx-1, Nx-1), speye(Nx-1));
  F = BuildCNRHS(Nx, ht, Nt, yhat(x(2:end-1), t', t(end), nu), y0(x(2:end-1)), sparse(Nx-1, 1));
  U = M\F; U = reshape(U, Nx-1, 2*Nt+2);
  y = U(:, 1:Nt+1); p = U(:, Nt+2:end);

  %same space time matrix for each subdomain
  Mm = BuildCNMatrix(A, Nx, ht, Twl/ht, nu, speye(Nx-1), sparse(Nx-1, Nx-1), sparse(Nx-1, Nx-1), speye(Nx-1)); 
  dMm = decomposition(Mm, 'lu');

  %transmission conditions (tc)
  gl = zeros(Nx-1, Ntw(k)-1); gr = gl; 

  for nit = 1:Nit
    %temporal tc and sol
    gltp = zeros(size(gl)); grtp = gltp;
    ytp = [];ptp = [];

    for m = 1:Ntw(k)
      if m == 1
        %build rhs and solve in Q1
        F1 = BuildCNRHS(Nx, ht, Twl/ht, yhat(x(2:end-1), t(1:Twl/ht+1)', t(end), nu), y0(x(2:end-1)), gr(:, 1));               
        U1 = Mm\F1; U1 = reshape(U1, Nx-1, 2*Twl/ht+2);
        y1 = U1(:, 1:Twl/ht+1); p1 = U1(:, Twl/ht+2:end);
            
        %assemble solutions and tc in Q1
        ytp = [ytp y1]; ptp = [ptp p1(:, 1:end-1)];
        gltp(:, 1) = y1(:, end);

      elseif m == Ntw(k)
        %build rhs and solve in Qn
        Fn = BuildCNRHS(Nx, ht, Twl/ht, yhat(x(2:end-1), t((Ntw(k)-1)*Twl/ht+1:end)', t(end), nu), gl(:, end), sparse(Nx-1, 1));        
        Un = Mm\Fn; Un = reshape(Un, Nx-1, 2*Twl/ht+2);
        yn = Un(:, 1:Twl/ht+1); pn = Un(:, Twl/ht+2:end);

        %assemble solutions and tc in Qn
        ytp = [ytp yn(:, 2:end)]; ptp = [ptp pn];
        grtp(:, end) = pn(:, 1);

      else
        %build rhs and solve in Qm               
        Fm = BuildCNRHS(Nx, ht, Twl/ht, yhat(x(2:end-1), t((m-1)*Twl/ht+1:m*Twl/ht+1)', t(end), nu), gl(:, m-1), gr(:, m));
        Um = Mm\Fm; Um = reshape(Um, Nx-1, 2*Twl/ht+2);
        ym = Um(:, 1:Twl/ht+1); pm = Um(:, Twl/ht+2:end);

        %assemble solutions and tc in Qm
        ytp = [ytp ym(:, 2:end)]; ptp = [ptp pm(:, 1:end-1)];
        gltp(:, m) = ym(:, end); grtp(:, m-1) = pm(:, 1);   
      end
    end

    %update tc and compute L^2 error
    gl = gltp; gr = grtp;                   
    erry(k,nit) = norm(ht*(full(y)-ytp));
    errp(k,nit) = norm(ht*(full(p)-ptp));   
  end
end

markers = {'*','s','+','o','x','d','p','v','h'};
nMarkers = length(markers);

figure
for k = 1:length(Ntw)
  semilogy(1:Nit, (erry(k, :)+errp(k, :))/(erry(k, 1)+errp(k, 1)), '--', 'Marker', markers{mod(k-1,length(markers))+1}, 'DisplayName', sprintf('$$N = $$ %d', Ntw(k)), 'linewidth', 2, 'MarkerSize', 15);
  hold on
end
semilogy(1:Nit, boundT(1).^(1:Nit), 'DisplayName', sprintf('Theo. Pred.'), 'linewidth', 2);
hold off
set(gca, 'FontSize', 15); 
set(gca, 'linewidth', 1.5);
grid on
ax = gca;
ax.YMinorGrid = 'off';
xlim([1 Nit]); ylim([1e-10 1e0]);
xlabel('Iteration', 'interpreter', 'latex', 'FontSize', 20);
ylabel('$$L^2$$ error', 'interpreter', 'latex', 'FontSize', 20);
legend('interpreter', 'latex', 'location', 'southwest', 'FontSize', 20);
title(['Time window length = 1/', num2str(1/Twl)], 'interpreter', 'latex', 'FontSize', 20);

%saveas(gca, 'HeatWS1d_nu_1e-1_L_1_hxht_2^-5_Twl_2^-1_Nit_10_Ntw_2-16', 'epsc');


figure
plot(2:Ntw(end), rhoT, 'k-o', 'linewidth', 2, 'Markersize', 15); 
hold on
plot(2:Ntw(end), normT, 'k--s', 'linewidth', 2, 'Markersize', 15); 
plot(2:Ntw(end), boundT, 'k-.^', 'linewidth', 2, 'Markersize', 15); 
hold off
set(gca, 'FontSize', 15); 
set(gca, 'linewidth', 1.5);
grid on
xlabel('$N$', 'interpreter', 'latex', 'FontSize', 20);
xlim([2 Ntw(end)]);
legend({'$$\rho(T^{PS}_{N,m})$$', '$$\|T^{PS}_{N,m}\|_{\infty}$$', '$$\widetilde \rho(m)$$'}, 'interpreter', 'latex', 'FontSize', 20, 'location', 'best');
title(['Time window length = 1/', num2str(1/Twl)], 'interpreter', 'latex', 'FontSize', 20);

%saveas(gca, 'CompareThreeCurves_nu_1e-1_L_1_hxht_2^-5_Twl_2^-1_Ntw_2-16', 'epsc');