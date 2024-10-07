function [delay_1, delay_3, total_delay, energy_1, energy_3, energy_initial] = systemupdate_fg(f, g, y1, y2, delay_2, energy_2, kappa)
%SYSTEMUPDATE_FG 此处显示有关此函数的摘要
%   此处显示详细说明

delay_1 = y1 ./ f;
delay_3 = y2 ./ g;
total_delay = delay_1 + delay_2 + delay_3;

energy_1 = kappa * y1 .* f.^2;
energy_3 = kappa * y2 .* g.^2;
total_energy = energy_1 + energy_2 + energy_3;
energy_initial = sum(total_energy);
end

