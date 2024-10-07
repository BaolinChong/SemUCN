function [gamma_k, Rate_k, delay_2, total_delay, energy_2, energy_initial] = systemupdate_pd(rho_k, p_d, GLk, sigma_n, N, Bandwidth, lambda, beta, eta, delay_1, delay_3, energy_1, energy_3)
%SYSTEMUPDATE_PD 此处显示有关此函数的摘要
%   此处显示详细说明

gamma_k = N * (sum(sqrt(p_d .* lambda))).^2 ./ ( sum(p_d' * beta) + 1);
Rate_k = Bandwidth * (1-eta) * log2(1 + gamma_k');
delay_2 = GLk .* rho_k ./ Rate_k;
total_delay = delay_1 + delay_2 + delay_3;

energy_2 = delay_2 .* (sum(p_d)' * sigma_n);
total_energy = energy_1 + energy_2 + energy_3;
energy_initial = sum(total_energy);

end

