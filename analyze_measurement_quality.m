function analyze_measurement_quality()
% This function analyzes the impact of measurement quality (std dev) on the
% estimation performance for Doppler-only, AoA-only, and combined approaches.


%fixed parameters
g = 3;
alt_kft = 10;
vel = 200;
p_true = [5000, 5000, 0, 1e9];
p_est_init = [5100, 4900, 0, 1.000001e9];
fo = 1e9;
L = 2;
T = 60;
del_T = 1;

% generate platform trajectory once
[Px, Py, Pz, Vx, Vy, Vz, mu_vect] = weave(g, T, del_T, alt_kft, vel);
t = 0:del_T:T;
Plat_Nav_Data = [Px; Py; Pz; Vx; Vy; Vz; t];

% different measurement quality parameters
dop_std_dev_values = logspace(-1, 1, 5);  % 0.1 to 10 Hz
aoa_std_dev_values = logspace(-1, 1, 5);  % 0.1 to 10 degrees

% store results
crlb_doppler = zeros(length(dop_std_dev_values), 1);
crlb_aoa = zeros(length(aoa_std_dev_values), 1);
crlb_combined = zeros(length(dop_std_dev_values), length(aoa_std_dev_values));

% prepare AoA-specific true and init values
p_true_aoa = p_true;
p_true_aoa(4) = 0;  % phase offset for AoA
p_est_init_aoa = p_est_init;
p_est_init_aoa(4) = 0;  % initial phase guess

% doppler only 
for i = 1:length(dop_std_dev_values)
    dop_std_dev = dop_std_dev_values(i);
    [~, H_true_dop] = doppler_loc(dop_std_dev, Plat_Nav_Data, p_true, p_est_init);
    J_dop = (1/(dop_std_dev^2)) * (H_true_dop' * H_true_dop);
    CRLB_dop = inv(J_dop);
    crlb_doppler(i) = trace(CRLB_dop(1:3, 1:3));
end

% AoA only
for j = 1:length(aoa_std_dev_values)
    aoa_std_dev = aoa_std_dev_values(j);
    aoa_std_dev_rad = aoa_std_dev * pi/180;
    [~, H_true_aoa] = LBI_loc(aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true_aoa, p_est_init_aoa, fo, L);
    J_aoa = (1/(aoa_std_dev_rad^2)) * (H_true_aoa' * H_true_aoa);
    CRLB_aoa = inv(J_aoa);
    crlb_aoa(j) = trace(CRLB_aoa(1:3, 1:3));
end

% combined approach - evaluate all combinations
for i = 1:length(dop_std_dev_values)
    dop_std_dev = dop_std_dev_values(i);
    for j = 1:length(aoa_std_dev_values)
        aoa_std_dev = aoa_std_dev_values(j);
        aoa_std_dev_rad = aoa_std_dev * pi/180;
        
        [~, H_true_combined] = doppler_aoa_loc(dop_std_dev, aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L);
        C_dop = (dop_std_dev^2) * eye(length(Px));
        C_aoa = (aoa_std_dev_rad^2) * eye(length(Px));
        C_combined = blkdiag(C_dop, C_aoa);
        J_combined = H_true_combined' * inv(C_combined) * H_true_combined;
        CRLB_combined = inv(J_combined);
        crlb_combined(i, j) = trace(CRLB_combined(1:3, 1:3));
    end
end

% plot results
figure;

% doppler-only
subplot(3, 1, 1);
loglog(dop_std_dev_values, crlb_doppler, 'bo-', 'LineWidth', 2);
title('Doppler-only Estimation Performance vs. Measurement Quality');
xlabel('Doppler Measurement Std. Dev. (Hz)');
ylabel('CRLB Trace (log scale)');
grid on;
set(gcf, 'Renderer', 'painters');
drawnow;

% AoA-only
subplot(3, 1, 2);
loglog(aoa_std_dev_values, crlb_aoa, 'go-', 'LineWidth', 2);
title('AoA-only Estimation Performance vs. Measurement Quality');
xlabel('AoA Measurement Std. Dev. (degrees)');
ylabel('CRLB Trace (log scale)');
grid on;
set(gcf, 'Renderer', 'painters');
drawnow;

% combined (surface plot)
subplot(3, 1, 3);
[X, Y] = meshgrid(aoa_std_dev_values, dop_std_dev_values);
surf(X, Y, crlb_combined);
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log');
title('Combined Doppler-AoA Estimation Performance');
xlabel('AoA Std. Dev. (degrees)');
ylabel('Doppler Std. Dev. (Hz)');
zlabel('CRLB Trace (log scale)');
colorbar;
view(45, 30);
set(gcf, 'Renderer', 'painters');
drawnow;

% calculate improvement ratios
figure;
improvement_matrix = zeros(length(dop_std_dev_values), length(aoa_std_dev_values));
for i = 1:length(dop_std_dev_values)
    for j = 1:length(aoa_std_dev_values)
        % use the minimum of Doppler-only and AoA-only as the baseline
        baseline = min(crlb_doppler(i), crlb_aoa(j));
        improvement_matrix(i, j) = baseline / crlb_combined(i, j);
    end
end

% plot improvement ratio
surf(X, Y, improvement_matrix);
set(gca, 'XScale', 'log', 'YScale', 'log');
title('Improvement Ratio of Combined Approach over Best Individual Method');
xlabel('AoA Std. Dev. (degrees)');
ylabel('Doppler Std. Dev. (Hz)');
zlabel('Improvement Ratio');
colorbar;
view(45, 30);
set(gcf, 'Renderer', 'painters');
drawnow;

% display summary
fprintf('Measurement Quality Analysis Summary:\n');
fprintf('-----------------------------------\n');
[max_improvement, idx] = max(improvement_matrix(:));
[i_max, j_max] = ind2sub(size(improvement_matrix), idx);
fprintf('Best improvement occurs at:\n');
fprintf('  - Doppler Std Dev = %.2f Hz\n', dop_std_dev_values(i_max));
fprintf('  - AoA Std Dev = %.2f degrees\n', aoa_std_dev_values(j_max));
fprintf('Maximum improvement ratio: %.2f times\n', max_improvement);

% create a table of improvement ratios
fprintf('\nImprovement Ratio Table (rows: Doppler std dev, columns: AoA std dev):\n');
fprintf('Doppler/AoA |');
for j = 1:length(aoa_std_dev_values)
    fprintf(' %.2f deg  |', aoa_std_dev_values(j));
end
fprintf('\n');
fprintf('------------|');
for j = 1:length(aoa_std_dev_values)
    fprintf('-----------|');
end
fprintf('\n');

for i = 1:length(dop_std_dev_values)
    fprintf('%.2f Hz     |', dop_std_dev_values(i));
    for j = 1:length(aoa_std_dev_values)
        fprintf(' %.2f times |', improvement_matrix(i, j));
    end
    fprintf('\n');
end
end
