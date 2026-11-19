% --------------------------------------------- %
% Script to generate the figure of the weak
% scalability test for a 1D advection-diffusion 
% OCP. If the data is not available, run the
% script AdvDiff1dWS.m to first generate data.
% --------------------------------------------- %
clc; clear; close all

if isfile('AdvDiff1dWSdata2-512_nu_1e-1_c_1_kap_1e-2_L_1_hxht_2^-5_Twl_2^-1_Nit_20.mat') == 0
  disp('Data does not exists!');
  return;
end

load('AdvDiff1dWSdata2-512_nu_1e-1_c_1_kap_1e-2_L_1_hxht_2^-5_Twl_2^-1_Nit_20.mat');


markers = {'*','s','+','o','x','d','p','v','h'};
nMarkers = length(markers);

figure
for k = 1:length(Ntw)
    semilogy(1:Nit, (erry(k, :)+errp(k, :))/(erry(k, 1)+errp(k, 1)), '--', 'Marker', markers{mod(k-1,length(markers))+1}, 'DisplayName', sprintf('$N = $ %d', Ntw(k)), 'linewidth', 2, 'MarkerSize', 15);
    hold on
end
hold off
set(gca, 'FontSize', 15); set(gca, 'linewidth', 1.5);
xlim([1 Nit]);
grid on
xlabel('Iteration', 'interpreter', 'latex', 'FontSize', 20);
ylabel('$L^2$ error', 'interpreter', 'latex', 'FontSize', 20);
legend('interpreter', 'latex', 'location', 'southwest', 'FontSize', 20);

%saveas(gca, 'AdvDiff1dWS_nu_1e-1_c_1_kap_1e-2_L_1_hxht_2^-5_Twl_2^-1_Nit_20_2-512', 'epsc');