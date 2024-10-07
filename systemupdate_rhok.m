function [Sk, y1, y2, delay_1, delay_2, delay_3, total_delay, energy_1, energy_2, energy_3, energy_initial] = systemupdate_rhok(rho_k, f, g, p_d, Rate_k, Lk, GLk, kappa, sigma_n, C1, C2, C3, C4, C5)
%SYSTEMUPDATE_RHOK 此处显示有关此函数的摘要
%   此处显示详细说明
Sk = GLk .* rho_k; % Extraction data
y1 = compute_y1(Lk, rho_k, C1, C2, C3);
y2 = compute_y2(rho_k, C4, C5);
delay_1 = y1 ./ f;
delay_2 = GLk .* rho_k ./ Rate_k;
delay_3 = y2 ./ g;
total_delay = delay_1 + delay_2 + delay_3;

energy_1 = kappa * y1 .* f.^2;
energy_2 = delay_2 .* (sum(p_d)' * sigma_n);
energy_3 = kappa * y2 .* g.^2;
total_energy = energy_1 + energy_2 + energy_3;
energy_initial = sum(total_energy);
end

