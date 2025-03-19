function CRLB = compute_CRLB(method, H, std_dev)
% Computes the Cramer-Rao Lower Bound (CRLB) for different methods.
%
% Usage: CRLB = compute_CRLB(method, H, std_dev);
%
% Inputs:   method = 'doppler', 'aoa', or 'combined'
%           H = Jacobian matrix for the respective method
%           std_dev = standard deviation of measurements
%                     For 'combined', std_dev should be [std_dev_dop, std_dev_aoa]
%
% Outputs:  CRLB = Cramer-Rao Lower Bound matrix

switch lower(method)
    case 'doppler'
        % Compute FIM and CRLB for Doppler-only
        J = (1/(std_dev^2)) * (H' * H);
        CRLB = inv(J);
        
    case 'aoa'
        % Compute FIM and CRLB for AoA-only
        J = (1/(std_dev^2)) * (H' * H);
        CRLB = inv(J);
        
    case 'combined'
        % For combined method, std_dev is [std_dev_dop, std_dev_aoa]
        % H is [H_dop; H_aoa]
        std_dev_dop = std_dev(1);
        std_dev_aoa = std_dev(2);
        
        % Extract Doppler and AoA parts of H
        num_dop_meas = size(H, 1) / 2;
        H_dop = H(1:num_dop_meas, :);
        H_aoa = H(num_dop_meas+1:end, :);
        
        % Create covariance matrices
        C_dop = (std_dev_dop^2) * eye(num_dop_meas);
        C_aoa = (std_dev_aoa^2) * eye(num_dop_meas);
        C_combined = blkdiag(C_dop, C_aoa);
        
        % Compute FIM and CRLB
        J = H' * inv(C_combined) * H;
        CRLB = inv(J);
        
    otherwise
        error('Unknown method. Use "doppler", "aoa", or "combined".');
end

end