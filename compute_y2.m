function y2 = compute_y2(rho_k, C4, C5)
%COMPUTE_Y3 此处显示有关此函数的摘要
%   此处显示详细说明
% C4 = 1e6;
% C5 = 1.3;
y2 = C4 .* rho_k .^ (-C5);
end

