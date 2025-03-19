function run_doppler_aoa(g, T, del_T, dop_std_dev, aoa_std_dev, alt_kft, vel, p_true, p_est_init, fo, L)
%
% Usage: run_doppler_aoa(g, T, del_T, dop_std_dev, aoa_std_dev, alt_kft, vel, p_true, p_est_init, fo, L);
% 
% Inputs:   g = platform turn acceleration in g's (e.g., g = 2 is a 2g turn)
%           T = total time of observation in seconds
%           del_T = measurement time spacing in seconds
%           dop_std_dev = standard deviation of frequency estimates in Hz(scalar)
%           aoa_std_dev = standard deviation of AoA measurements in degrees
%           alt_kft = altitude of the platform in kft
%           vel = speed of platform in m/s
%           p_true = True Emitter Location & Frequency (1x4 Vector)
%                   [xe ye ze fo] with units [m m m Hz]
%           p_est_init = Initial Estimate of Emitter (1x4 Vector)
%                   [xe_init ye_init ze_init fo_init] with units [m m m Hz]
%           fo = signal's center frequency in Hz
%           L = baseline length in meters
%
% Outputs:  None - Creates figures showing results

% Generate Nav Data:
[Px, Py, Pz, Vx, Vy, Vz, mu_vect] = weave(g, T, del_T, alt_kft, vel);
t = 0:del_T:T;
Plat_Nav_Data = [Px; Py; Pz; Vx; Vy; Vz; t];

% Convert AoA standard deviation from degrees to radians
aoa_std_dev_rad = aoa_std_dev * pi/180; 

% run Monte Carlo Sims of Combined Doppler-AoA Location:
for n = 1:200
    [p_est_iter, H_true_combined] = doppler_aoa_loc(dop_std_dev, aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L);
    [~, N] = size(p_est_iter);
    loc_est(:, n) = p_est_iter(:, N);
end

% Calculate CRLB
C_dop = (dop_std_dev^2) * eye(length(Px));
C_aoa = (aoa_std_dev_rad^2) * eye(length(Px));
C_combined = blkdiag(C_dop, C_aoa);

J_combined = H_true_combined' * inv(C_combined) * H_true_combined;
CRLB_combined = inv(J_combined);
C_xy = CRLB_combined(1:2, 1:2);

% Calculate errors
x_error = loc_est(1, :) - p_true(1);
y_error = loc_est(2, :) - p_true(2);
z_error = loc_est(3, :) - p_true(3);

rms_x = sqrt(mean(x_error.^2));
rms_y = sqrt(mean(y_error.^2));
rms_z = sqrt(mean(z_error.^2));

% Plot results
figure;
plot(Px, Py, 'b');
hold on;
plot(loc_est(1, :), loc_est(2, :), 'o', p_true(1), p_true(2), 'rx');
title('Emitter Location Estimation using Combined Doppler and AoA');
xlabel('X Position (m)');
ylabel('Y Position (m)');
set(gcf, 'Renderer', 'painters');
drawnow;

hold on;
ellipse(C_xy, 6, 'r', p_true(1:2)); % using k=6 here gives 95% ellipse
hold off;

% The following 4 lines duplicates figure 1 into figure 2
h1 = gcf;
h2 = figure;
objects = allchild(h1);
copyobj(get(h1, 'children'), h2);

% Zoom in to show close up of ellipse
xe = p_true(1);
ye = p_true(2);
axis([[-6 6]*sqrt(C_xy(1,1))+xe  [-6 6]*sqrt(C_xy(2,2))+ye]) 
axis equal;
title('Zoomed View of Error Ellipse (95% Confidence)');
set(gcf, 'Renderer', 'painters');
drawnow;

% Display RMS errors
fprintf('\nCombined Doppler-AoA Results:\n');
fprintf('-------------------------\n');
fprintf('RMS Error in X: %.2f m\n', rms_x);
fprintf('RMS Error in Y: %.2f m\n', rms_y);
fprintf('RMS Error in Z: %.2f m\n', rms_z);
fprintf('3D RMS Error: %.2f m\n', sqrt(rms_x^2 + rms_y^2 + rms_z^2));
fprintf('CRLB Trace (Position): %.2e\n', trace(CRLB_combined(1:3, 1:3)));

end
