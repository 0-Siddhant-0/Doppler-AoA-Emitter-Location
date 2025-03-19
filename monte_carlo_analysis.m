function monte_carlo_analysis(g, T, del_T, dop_std_dev, aoa_std_dev, alt_kft, vel, p_true, p_est_init, fo, L, num_runs)
% This function performs a detailed Monte Carlo analysis to evaluate the 
% statistical performance of the three estimation methods (Doppler-only, 
% AoA-only, and combined).
%
% Usage: monte_carlo_analysis(g, T, del_T, dop_std_dev, aoa_std_dev, alt_kft, vel, p_true, p_est_init, fo, L, num_runs);
% 
% Inputs:   g = platform turn acceleration in g's (e.g., g = 2 is a 2g turn)
%           T = total time of observation in seconds
%           del_T = measurement time spacing in seconds
%           dop_std_dev = standard deviation of frequency estimates in Hz
%           aoa_std_dev = standard deviation of AoA measurements in degrees
%           alt_kft = altitude of the platform in kft
%           vel = speed of platform in m/s
%           p_true = True Emitter Location & Frequency (1x4 Vector)
%                   [xe ye ze fo] with units [m m m Hz]
%           p_est_init = Initial Estimate of Emitter (1x4 Vector)
%                   [xe_init ye_init ze_init fo_init] with units [m m m Hz]
%           fo = signal's center frequency in Hz
%           L = baseline length in meters
%           num_runs = number of Monte Carlo runs to perform
%
% Outputs: None (creates plots and prints statistics)

% If num_runs is not provided, use default value
if nargin < 12
    num_runs = 1000;  % Larger number for more statistically significant results
end

% Convert AoA standard deviation from degrees to radians
aoa_std_dev_rad = aoa_std_dev * pi/180;

% Generate platform trajectory
[Px, Py, Pz, Vx, Vy, Vz, mu_vect] = weave(g, T, del_T, alt_kft, vel);
t = 0:del_T:T;
Plat_Nav_Data = [Px; Py; Pz; Vx; Vy; Vz; t];

% Setup for AoA method
p_true_aoa = p_true;
p_true_aoa(4) = 0;  % phase offset for AoA
p_est_init_aoa = p_est_init;
p_est_init_aoa(4) = 0;  % initial phase guess

% Pre-allocate arrays for results
loc_est_dop = zeros(4, num_runs);
loc_est_aoa = zeros(4, num_runs);
loc_est_combined = zeros(4, num_runs);

% Progress reporting
fprintf('Running Monte Carlo Simulation with %d runs...\n', num_runs);
progress_step = max(1, floor(num_runs/20));  % Report progress in 5% increments

% Run Monte Carlo simulations
for n = 1:num_runs
    % Progress reporting
    if mod(n, progress_step) == 0
        fprintf('  Progress: %.1f%% complete (%d/%d runs)\n', 100*n/num_runs, n, num_runs);
    end
    
    % Doppler-only
    [p_est_iter_dop, ~] = doppler_loc(dop_std_dev, Plat_Nav_Data, p_true, p_est_init);
    [~, N_dop] = size(p_est_iter_dop);
    loc_est_dop(:, n) = p_est_iter_dop(:, N_dop);
    
    % AoA-only
    [p_est_iter_aoa, ~] = LBI_loc(aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true_aoa, p_est_init_aoa, fo, L);
    [~, N_aoa] = size(p_est_iter_aoa);
    loc_est_aoa(:, n) = p_est_iter_aoa(:, N_aoa);
    
    % Combined
    [p_est_iter_combined, ~] = doppler_aoa_loc(dop_std_dev, aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L);
    [~, N_combined] = size(p_est_iter_combined);
    loc_est_combined(:, n) = p_est_iter_combined(:, N_combined);
end

fprintf('  Progress: 100.0%% complete (%d/%d runs)\n', num_runs, num_runs);

% Compute estimation errors
x_err_dop = loc_est_dop(1, :) - p_true(1);
y_err_dop = loc_est_dop(2, :) - p_true(2);
z_err_dop = loc_est_dop(3, :) - p_true(3);
fo_err_dop = loc_est_dop(4, :) - p_true(4);

x_err_aoa = loc_est_aoa(1, :) - p_true(1);
y_err_aoa = loc_est_aoa(2, :) - p_true(2);
z_err_aoa = loc_est_aoa(3, :) - p_true(3);
fo_err_aoa = loc_est_aoa(4, :) - p_true(4);

x_err_combined = loc_est_combined(1, :) - p_true(1);
y_err_combined = loc_est_combined(2, :) - p_true(2);
z_err_combined = loc_est_combined(3, :) - p_true(3);
fo_err_combined = loc_est_combined(4, :) - p_true(4);

% Calculate RMS errors
rms_x_dop = sqrt(mean(x_err_dop.^2));
rms_y_dop = sqrt(mean(y_err_dop.^2));
rms_z_dop = sqrt(mean(z_err_dop.^2));
rms_fo_dop = sqrt(mean(fo_err_dop.^2));

rms_x_aoa = sqrt(mean(x_err_aoa.^2));
rms_y_aoa = sqrt(mean(y_err_aoa.^2));
rms_z_aoa = sqrt(mean(z_err_aoa.^2));
rms_fo_aoa = sqrt(mean(fo_err_aoa.^2));

rms_x_combined = sqrt(mean(x_err_combined.^2));
rms_y_combined = sqrt(mean(y_err_combined.^2));
rms_z_combined = sqrt(mean(z_err_combined.^2));
rms_fo_combined = sqrt(mean(fo_err_combined.^2));

% Compute 3D position errors
pos_err_dop = sqrt(x_err_dop.^2 + y_err_dop.^2 + z_err_dop.^2);
pos_err_aoa = sqrt(x_err_aoa.^2 + y_err_aoa.^2 + z_err_aoa.^2);
pos_err_combined = sqrt(x_err_combined.^2 + y_err_combined.^2 + z_err_combined.^2);

% Calculate CRLB for reference
[~, H_true_dop] = doppler_loc(dop_std_dev, Plat_Nav_Data, p_true, p_est_init);
J_dop = (1/(dop_std_dev^2)) * (H_true_dop' * H_true_dop);
CRLB_dop = inv(J_dop);

[~, H_true_aoa] = LBI_loc(aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true_aoa, p_est_init_aoa, fo, L);
J_aoa = (1/(aoa_std_dev_rad^2)) * (H_true_aoa' * H_true_aoa);
CRLB_aoa = inv(J_aoa);

[~, H_true_combined] = doppler_aoa_loc(dop_std_dev, aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L);
C_dop = (dop_std_dev^2) * eye(length(Px));
C_aoa = (aoa_std_dev_rad^2) * eye(length(Px));
C_combined = blkdiag(C_dop, C_aoa);
J_combined = H_true_combined' * inv(C_combined) * H_true_combined;
CRLB_combined = inv(J_combined);

% Calculate the theoretical position RMS errors based on CRLB
crlb_rms_pos_dop = sqrt(trace(CRLB_dop(1:3, 1:3)));
crlb_rms_pos_aoa = sqrt(trace(CRLB_aoa(1:3, 1:3)));
crlb_rms_pos_combined = sqrt(trace(CRLB_combined(1:3, 1:3)));

% Plot position error histograms
figure;
edges = linspace(0, max([max(pos_err_dop), max(pos_err_aoa), max(pos_err_combined)]), 30);

subplot(3, 1, 1);
histogram(pos_err_dop, edges, 'FaceColor', 'b');
hold on;
xline(sqrt(mean(pos_err_dop.^2)), 'r', 'LineWidth', 2);
xline(crlb_rms_pos_dop, 'g--', 'LineWidth', 2);
title('Doppler-only Position Error Distribution');
xlabel('3D Position Error (m)');
ylabel('Frequency');
legend('Error Distribution', 'RMS Error', 'CRLB Bound');
grid on;
set(gca, 'Renderer', 'painters');
drawnow;

subplot(3, 1, 2);
histogram(pos_err_aoa, edges, 'FaceColor', 'g');
hold on;
xline(sqrt(mean(pos_err_aoa.^2)), 'r', 'LineWidth', 2);
xline(crlb_rms_pos_aoa, 'g--', 'LineWidth', 2);
title('AoA-only Position Error Distribution');
xlabel('3D Position Error (m)');
ylabel('Frequency');
legend('Error Distribution', 'RMS Error', 'CRLB Bound');
grid on;
set(gca, 'Renderer', 'painters');
drawnow;

subplot(3, 1, 3);
histogram(pos_err_combined, edges, 'FaceColor', 'r');
hold on;
xline(sqrt(mean(pos_err_combined.^2)), 'r', 'LineWidth', 2);
xline(crlb_rms_pos_combined, 'g--', 'LineWidth', 2);
title('Combined Doppler-AoA Position Error Distribution');
xlabel('3D Position Error (m)');
ylabel('Frequency');
legend('Error Distribution', 'RMS Error', 'CRLB Bound');
grid on;
set(gca, 'Renderer', 'painters');
drawnow;

% Plot scatter of position estimates
figure;
plot3(loc_est_dop(1, :), loc_est_dop(2, :), loc_est_dop(3, :), 'b.', 'MarkerSize', 3);
hold on;
plot3(loc_est_aoa(1, :), loc_est_aoa(2, :), loc_est_aoa(3, :), 'g.', 'MarkerSize', 3);
plot3(loc_est_combined(1, :), loc_est_combined(2, :), loc_est_combined(3, :), 'r.', 'MarkerSize', 3);
plot3(p_true(1), p_true(2), p_true(3), 'k*', 'MarkerSize', 10, 'LineWidth', 2);
title('3D Scatter Plot of Location Estimates');
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
legend('Doppler-only', 'AoA-only', 'Combined', 'True Location');
grid on;
view(45, 30);
set(gca, 'Renderer', 'painters');
drawnow;

% Print detailed statistical summary
fprintf('\nMonte Carlo Analysis Results (%d runs):\n', num_runs);
fprintf('-------------------------------------\n');
fprintf('\nRMS Position Errors:\n');
fprintf('Method     | X Error (m) | Y Error (m) | Z Error (m) | 3D Error (m) | Freq Error (Hz)\n');
fprintf('Doppler    | %.2f        | %.2f        | %.2f        | %.2f         | %.2f\n', ...
    rms_x_dop, rms_y_dop, rms_z_dop, sqrt(rms_x_dop^2 + rms_y_dop^2 + rms_z_dop^2), rms_fo_dop);
fprintf('AoA        | %.2f        | %.2f        | %.2f        | %.2f         | %.2f\n', ...
    rms_x_aoa, rms_y_aoa, rms_z_aoa, sqrt(rms_x_aoa^2 + rms_y_aoa^2 + rms_z_aoa^2), rms_fo_aoa);
fprintf('Combined   | %.2f        | %.2f        | %.2f        | %.2f         | %.2f\n', ...
    rms_x_combined, rms_y_combined, rms_z_combined, sqrt(rms_x_combined^2 + rms_y_combined^2 + rms_z_combined^2), rms_fo_combined);

fprintf('\nTheoretical CRLB Position Errors:\n');
fprintf('Method     | CRLB RMS (m) | Actual/CRLB Ratio\n');
fprintf('Doppler    | %.2f         | %.2f\n', ...
    crlb_rms_pos_dop, sqrt(mean(pos_err_dop.^2))/crlb_rms_pos_dop);
fprintf('AoA        | %.2f         | %.2f\n', ...
    crlb_rms_pos_aoa, sqrt(mean(pos_err_aoa.^2))/crlb_rms_pos_aoa);
fprintf('Combined   | %.2f         | %.2f\n', ...
    crlb_rms_pos_combined, sqrt(mean(pos_err_combined.^2))/crlb_rms_pos_combined);

fprintf('\nPerformance Improvement Ratios:\n');
fprintf('Doppler -> Combined: %.2f times improvement\n', ...
    sqrt(mean(pos_err_dop.^2))/sqrt(mean(pos_err_combined.^2)));
fprintf('AoA -> Combined: %.2f times improvement\n', ...
    sqrt(mean(pos_err_aoa.^2))/sqrt(mean(pos_err_combined.^2)));

% Calculate circular error probable (CEP) statistics
pos_err_dop_sorted = sort(pos_err_dop);
pos_err_aoa_sorted = sort(pos_err_aoa);
pos_err_combined_sorted = sort(pos_err_combined);

cep50_dop = pos_err_dop_sorted(ceil(0.5*num_runs));
cep90_dop = pos_err_dop_sorted(ceil(0.9*num_runs));
cep95_dop = pos_err_dop_sorted(ceil(0.95*num_runs));

cep50_aoa = pos_err_aoa_sorted(ceil(0.5*num_runs));
cep90_aoa = pos_err_aoa_sorted(ceil(0.9*num_runs));
cep95_aoa = pos_err_aoa_sorted(ceil(0.95*num_runs));

cep50_combined = pos_err_combined_sorted(ceil(0.5*num_runs));
cep90_combined = pos_err_combined_sorted(ceil(0.9*num_runs));
cep95_combined = pos_err_combined_sorted(ceil(0.95*num_runs));

fprintf('\nCircular Error Probable (CEP) Statistics:\n');
fprintf('Method     | CEP50 (m)   | CEP90 (m)   | CEP95 (m)\n');
fprintf('Doppler    | %.2f        | %.2f        | %.2f\n', cep50_dop, cep90_dop, cep95_dop);
fprintf('AoA        | %.2f        | %.2f        | %.2f\n', cep50_aoa, cep90_aoa, cep95_aoa);
fprintf('Combined   | %.2f        | %.2f        | %.2f\n', cep50_combined, cep90_combined, cep95_combined);

end
