function compare_methods(g, T, del_T, dop_std_dev, aoa_std_dev, alt_kft, vel, p_true, p_est_init, fo, L)
% This function performs a detailed comparison of Doppler-only, AoA-only, and 
% combined approaches using Monte Carlo simulations and CRLB analysis.
%
% Usage: compare_methods(g, T, del_T, dop_std_dev, aoa_std_dev, alt_kft, vel, p_true, p_est_init, fo, L);
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


% convert AoA standard deviation from degrees to radians
aoa_std_dev_rad = aoa_std_dev * pi/180;

% generate platform trajectory
[Px, Py, Pz, Vx, Vy, Vz, mu_vect] = weave(g, T, del_T, alt_kft, vel);
t = 0:del_T:T;
Plat_Nav_Data = [Px; Py; Pz; Vx; Vy; Vz; t];

% setup for AoA method
p_true_aoa = p_true;
p_true_aoa(4) = 0;  % phase offset for AoA
p_est_init_aoa = p_est_init;
p_est_init_aoa(4) = 0;  % initial phase guess

% number of Monte Carlo runs
num_mc = 10000;

% pre-allocate arrays
loc_est_dop = zeros(4, num_mc);
loc_est_aoa = zeros(4, num_mc);
loc_est_combined = zeros(4, num_mc);

% run Monte Carlo simulations
for n = 1:num_mc
    % Doppler-only
    [p_est_iter_dop, ~] = doppler_loc(dop_std_dev, Plat_Nav_Data, p_true, p_est_init);
    [~, N_dop] = size(p_est_iter_dop);
    loc_est_dop(:, n) = p_est_iter_dop(:, N_dop);
    
    % AoA-only
    [p_est_iter_aoa, ~] = LBI_loc(aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true_aoa, p_est_init_aoa, fo, L);
    [~, N_aoa] = size(p_est_iter_aoa);
    loc_est_aoa(:, n) = p_est_iter_aoa(:, N_aoa);
    
    % combined
    [p_est_iter_combined, ~] = doppler_aoa_loc(dop_std_dev, aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L);
    [~, N_combined] = size(p_est_iter_combined);
    loc_est_combined(:, n) = p_est_iter_combined(:, N_combined);
end

% compute CRLB for all methods
[~, H_true_dop] = doppler_loc(dop_std_dev, Plat_Nav_Data, p_true, p_est_init);
J_dop = (1/(dop_std_dev^2)) * (H_true_dop' * H_true_dop);
CRLB_dop = inv(J_dop);
C_xy_dop = CRLB_dop(1:2, 1:2);

[~, H_true_aoa] = LBI_loc(aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true_aoa, p_est_init_aoa, fo, L);
J_aoa = (1/(aoa_std_dev_rad^2)) * (H_true_aoa' * H_true_aoa);
CRLB_aoa = inv(J_aoa);
C_xy_aoa = CRLB_aoa(1:2, 1:2);

[~, H_true_combined] = doppler_aoa_loc(dop_std_dev, aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L);
C_dop = (dop_std_dev^2) * eye(length(Px));
C_aoa = (aoa_std_dev_rad^2) * eye(length(Px));
C_combined = blkdiag(C_dop, C_aoa);
J_combined = H_true_combined' * inv(C_combined) * H_true_combined;
CRLB_combined = inv(J_combined);
C_xy_combined = CRLB_combined(1:2, 1:2);

% calculate estimation errors
x_err_dop = loc_est_dop(1, :) - p_true(1);
y_err_dop = loc_est_dop(2, :) - p_true(2);
z_err_dop = loc_est_dop(3, :) - p_true(3);

x_err_aoa = loc_est_aoa(1, :) - p_true(1);
y_err_aoa = loc_est_aoa(2, :) - p_true(2);
z_err_aoa = loc_est_aoa(3, :) - p_true(3);

x_err_combined = loc_est_combined(1, :) - p_true(1);
y_err_combined = loc_est_combined(2, :) - p_true(2);
z_err_combined = loc_est_combined(3, :) - p_true(3);

% calculate RMS errors
rms_x_dop = sqrt(mean(x_err_dop.^2));
rms_y_dop = sqrt(mean(y_err_dop.^2));
rms_z_dop = sqrt(mean(z_err_dop.^2));

rms_x_aoa = sqrt(mean(x_err_aoa.^2));
rms_y_aoa = sqrt(mean(y_err_aoa.^2));
rms_z_aoa = sqrt(mean(z_err_aoa.^2));

rms_x_combined = sqrt(mean(x_err_combined.^2));
rms_y_combined = sqrt(mean(y_err_combined.^2));
rms_z_combined = sqrt(mean(z_err_combined.^2));

% create 2D scatter plot with error ellipses
figure;
hold on;
plot(p_true(1), p_true(2), 'k*', 'MarkerSize', 10);
plot(loc_est_dop(1, :), loc_est_dop(2, :), 'bx', 'MarkerSize', 3);
plot(loc_est_aoa(1, :), loc_est_aoa(2, :), 'gx', 'MarkerSize', 3);
plot(loc_est_combined(1, :), loc_est_combined(2, :), 'rx', 'MarkerSize', 3);
ellipse(C_xy_dop, 6, 'b', p_true(1:2));
ellipse(C_xy_aoa, 6, 'g', p_true(1:2));
ellipse(C_xy_combined, 6, 'r', p_true(1:2));
title('Comparison of Location Estimates and Error Ellipses (95% Confidence)');
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('True Location', 'Doppler-only', 'AoA-only', 'Combined Doppler-AoA');
set(gcf, 'Renderer', 'painters');
drawnow;
hold off;

% create duplicate figure for zoomed view
h1 = gcf;
h2 = figure;
copyobj(get(h1, 'children'), h2);

% zoomed view of error ellipses
title('Zoomed View of Error Ellipses (95% Confidence)');
% focus on the area around the true location
xe = p_true(1);
ye = p_true(2);
max_semi_major = max([sqrt(max(eig(C_xy_dop))), sqrt(max(eig(C_xy_aoa))), sqrt(max(eig(C_xy_combined)))]);
zoom_factor = 6;  % show 6x the largest semi-major axis
axis([xe-zoom_factor*max_semi_major, xe+zoom_factor*max_semi_major, ...
      ye-zoom_factor*max_semi_major, ye+zoom_factor*max_semi_major]);
axis equal;
set(gcf, 'Renderer', 'painters');
drawnow;

% print results summary
fprintf('\nMethods Comparison Summary:\n');
fprintf('---------------------------\n');
fprintf('\nRMS Position Errors:\n');
fprintf('Method     | X Error (m) | Y Error (m) | Z Error (m) | 3D Error (m)\n');
fprintf('Doppler    | %.2f        | %.2f        | %.2f        | %.2f\n', ...
    rms_x_dop, rms_y_dop, rms_z_dop, sqrt(rms_x_dop^2 + rms_y_dop^2 + rms_z_dop^2));
fprintf('AoA        | %.2f        | %.2f        | %.2f        | %.2f\n', ...
    rms_x_aoa, rms_y_aoa, rms_z_aoa, sqrt(rms_x_aoa^2 + rms_y_aoa^2 + rms_z_aoa^2));
fprintf('Combined   | %.2f        | %.2f        | %.2f        | %.2f\n', ...
    rms_x_combined, rms_y_combined, rms_z_combined, sqrt(rms_x_combined^2 + rms_y_combined^2 + rms_z_combined^2));

fprintf('\nCRLB Analysis:\n');
fprintf('Method     | CRLB Trace  | Improvement over Doppler | Improvement over AoA\n');
fprintf('Doppler    | %.2e   | -                        | -\n', trace(CRLB_dop(1:3, 1:3)));
fprintf('AoA        | %.2e   | -                        | -\n', trace(CRLB_aoa(1:3, 1:3)));
fprintf('Combined   | %.2e   | %.2f times                | %.2f times\n', ...
    trace(CRLB_combined(1:3, 1:3)), ...
    trace(CRLB_dop(1:3, 1:3))/trace(CRLB_combined(1:3, 1:3)), ...
    trace(CRLB_aoa(1:3, 1:3))/trace(CRLB_combined(1:3, 1:3)));

end
