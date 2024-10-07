function y1 = compute_y1(Lk, rho_k, C1, C2, C3)
%COMPUTE_Y1 此处显示有关此函数的摘要
%   D: Dk,bit
%   S: Sk,bit
% C1 = Lk * 5.5e1;   % 6.7
% C2 = 0.5;
% C3 = 2;
weight = 0.55;
% y1 = weight * Lk.^(1.25) + (1 - weight) * Lk.^(0.75) - C1 .* (rho_k - C2).^C3; % 5.5

y1 = weight * Lk.^(1.1) + (1 - weight) * Lk.^(0.9) - C1 .* (rho_k - C2).^C3 + max(C1 .* (0 - C2).^C3, C1 .* (1 - C2).^C3);

end

