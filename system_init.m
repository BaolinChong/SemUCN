% clc
% clear all
times_now = 1;
%% UCN network generation
% system parament intilization
M = 25; % the number of APs
% K = 16; % the number of UEs, the pilot length is equal to the number of UE
L = 200; % number of channel uses
eta = K/L;

% % variable
N = 36;
Delay_req = ones(K, 1) * 0.06; % performance requirement, s
rho_req = 0.1 * ones(K,1); % performance requirement, 0-1

rho_p = 0.1; % pilot signal power 20 dBw
rho_d = 0.2; % maximum downlink transmit power 23 dBw each AP

sigma_s = 8; % shadow fading dB
h_AP = 15; % height of AP, m
h_UE = 1.65; % height of user, m
f_carrier = 1.9e3; % carrier frequency, MHz
d_0 = 10; % m
d_1 = 50; % m
Bandwidth = 10e6; % system bandwidth in Hz
sigma_n = Bandwidth * 3.9810717055349565e-21; % noise variance -174dBm/Hz

rho_p = rho_p / sigma_n; % normlization pilot power
rho_d = rho_d / sigma_n; % normlization dowmlink transmit power

location_APs_load= cell2mat(struct2cell(load('./Location_information/location_information_rng10', 'location_APs')));
location_UEs_load= cell2mat(struct2cell(load('./Location_information/location_information_rng10', 'location_UEs')));
location_APs = reshape(location_APs_load(times_now,1:M,:),[M,2]); % 2-D location of APs we consider
location_UEs = reshape(location_UEs_load(times_now,1:K,:),[K,2]); % 2-D location of UEs we consider

% caculation the distance between APs and UE, M x K matrix
distance_matrix = zeros(M,K);
for m = 1:M
    for k = 1:K
        distance_matrix(m,k) = sqrt((location_APs(m,1) - location_UEs(k,1))^2 + (location_APs(m,2) - location_UEs(k,2))^2);
    end
end
% calculate large-scale fading
% shadow fading
a_shadow_load = cell2mat(struct2cell(load('./Location_information/location_information_rng10', 'a_shadow')));
b_shadow_load = cell2mat(struct2cell(load('./Location_information/location_information_rng10', 'b_shadow')));
a_shadow = a_shadow_load(times_now,1:M);
b_shadow = b_shadow_load(times_now,1:K);
delta_shadow_load = cell2mat(struct2cell(load('./Location_information/location_information_rng10', 'delta_shadow')));
delta_shadow = delta_shadow_load(times_now);
theta_shadow = zeros(M, K);
for m = 1:M
    for k = 1:K
        theta_shadow(m, k) = sigma_s * (sqrt(delta_shadow) * a_shadow(1, m) + sqrt(1 - delta_shadow) * b_shadow(1, k));
    end
end
% three-slope path loss 
L_loss = 46.3 + 33.9 * log10(f_carrier) - 13.82 * log10(h_AP) - (1.11 * log10(f_carrier) - 0.7) * h_UE + 1.56 * log10(f_carrier) - 0.8;
PL_loss = zeros(M, K);
for m = 1:M
    for k = 1:K
        if distance_matrix(m, k) > d_1
            PL_loss(m, k) = -L_loss - 35 * log10(distance_matrix(m, k)/1000) - theta_shadow(m, k);
        elseif (d_0 <= distance_matrix(m, k)) && (distance_matrix(m, k) <= d_1)
            PL_loss(m, k) = -L_loss - 10 * log10(d_1^1.5 * (distance_matrix(m, k)/1000)^2);
        else
            PL_loss(m, k) = -L_loss - 10 * log10(d_1^1.5 * d_0^2);
        end
    end
end
beta = 10.^(PL_loss / 10);

% uplink data extiamtion
lambda = K * rho_p * beta.^2 ./ (K * rho_p * beta + 1);

% intilization power allocation; large-scale-fading-based scheme
p_d = rho_d * beta ./ sum(beta,2);
% p_d_eval = [5023772863019.18] * 1/K * ones(M,K);

gamma_k = N * (sum(sqrt(p_d .* lambda))).^2 ./ ( sum(p_d' * beta) + 1);
Rate_k = Bandwidth * (1-eta) * log2(1 + gamma_k');

%% SemComm intilization
rng(10);
% Lk_max = 1100;    % from [500:100:1100]
Lk_min = 100;    % Kb
Lk = (Lk_min + (Lk_max - Lk_min) * rand(K, 1)) * 1000;  % b
weight_Lk = 0.5;
GLk = weight_Lk * Lk + (1 - weight_Lk) * Lk.^(0.8);

rho_k = 0.75 * ones(K,1); % Extraction rate intilization

f_max = 10e9;   % maximum compution capacity at CPU
g_max = 2e9;    % maximum compution capacity at User
f = ones(K, 1) * f_max / K;   % Intilization CPU capacity allocation
g = ones(K, 1) * g_max * 0.6;   % Intilization user capacity control

C1 = Lk * 3.0e1;  C2 = 0.65;  C3 = 2;
C4 = 1e6 + 2 * 1e5 * (2 * rand(K,1) - 1);  C5 = 1.3 + 2 * 0.13 * (2 * rand(K,1) - 1);
kappa = 2e-28; % power consumption coefficient

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



