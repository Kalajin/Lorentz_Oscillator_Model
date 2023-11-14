clear all

% Define the wavelength range
wavelength_um_range = 0.3:0.1:20; % Adjust the range as needed

%Define resonant wavenumber for each oscillator
wavenumber = [1180, 1108, 1063, 981, 634, 603 ]; % Input based on wavenumber cm-1

% Converts wavenumber to wavelength
resonant_wavelengths_um = (1./(wavenumber.*10.^-2)).*10.^2;

%Oscillator Strength scale
f0_scale=1;
gamma_scale=1e13;
omega_p_scale=3e13;
% Define oscillator parameters [f0, gamma, omega_p] for each oscillator
oscillator_params = [
    1*f0_scale, 0.5*gamma_scale, 0.9*omega_p_scale;
    1*f0_scale, 1.5*gamma_scale, 2.2*omega_p_scale;
    1*f0_scale, 0.8*gamma_scale, 1.5*omega_p_scale;
    1*f0_scale, 0.5*gamma_scale, 1*omega_p_scale;
    1*f0_scale, 0.45*gamma_scale, 0.8*omega_p_scale;
    1*f0_scale, 0.4*gamma_scale, 1*omega_p_scale;
];

% Calculate the complex refractive index
n_complex = lorentzOscillator(wavelength_um_range, resonant_wavelengths_um, oscillator_params);



% Plot the real and imaginary parts of the complex refractive index
figure;
subplot(2, 1, 1);
plot(wavelength_um_range, real(n_complex), 'b-');
xlabel('Wavelength (\mum)');
ylabel('Real Refractive Index (n)');
title('Refractive Index vs. Wavelength');

subplot(2, 1, 2);
plot(wavelength_um_range, imag(n_complex), 'r-');
xlabel('Wavelength (\mum)');
ylabel('Extinction Coefficient (k)');
title('Extinction Coefficient vs. Wavelength');

% Create a matrix with wavelength, n, and k values
results_matrix(:,1) = wavelength_um_range;
results_matrix(:,2) = real(n_complex);
results_matrix(:,3) = imag(n_complex);

% Save the results to an Excel file
filename = 'BaSO4 refractive_index_results_multi.xlsx';
writematrix(results_matrix, filename);

% Display a message indicating the file has been saved
fprintf('Results saved to %s\n', filename);


function n_complex = lorentzOscillator(wavelength_um_range, resonant_wavelengths_um, oscillator_params)
    % Input:
    % - wavelength_um_range: Wavelength range (in micrometers)
    % - resonant_wavelengths_um: Resonant wavelengths (in micrometers)
    % - oscillator_params: A matrix where each row contains [f0, gamma, omega_p] for each oscillator

    % Constants
    e_o = 8.854e-12;  % Permittivity of free space (F/m)
    c = 299792458;    % Speed of light (m/s)

    % Angular frequency
    wavelength_range_micrometers = wavelength_um_range;  % Ensure consistency
    omega = 2 * pi * c ./ (wavelength_range_micrometers * 1e-6);  % Convert to meters

    % Initialize complex refractive index
    n_complex = zeros(length(wavelength_range_micrometers), 1);

    % Calculate the complex refractive index using the Lorentz oscillator model
    for i = 1:length(wavelength_range_micrometers)
        w = omega(i);

        % Initialize real and imaginary parts
        n_real = 1.64;
        n_imag = 0;

        % Calculate contributions from each oscillator
        for m = 1:size(oscillator_params, 1)
            f0 = oscillator_params(m, 1);
            gamma = oscillator_params(m, 2);
            omega_p = oscillator_params(m, 3);
            w0 = 2 * pi * c ./ (resonant_wavelengths_um(m) * 1e-6);  % Convert to meters

            % Check for physically meaningful parameters
%             if w0 > w
%                 continue;  % Skip this oscillator if the resonant frequency is higher than the angular frequency
%             end

            % Contributions to n_real and n_imag for each oscillator
            n_real = n_real + (f0 * omega_p^2 * (w0^2 - w^2)) / ((w0^2 - w^2)^2 + w^2 * gamma^2);
            n_imag = n_imag + (f0 * omega_p^2 * gamma * w) / ((w0^2 - w^2)^2 + w^2 * gamma^2);
        end

        % Complex refractive index at this wavelength
        n_complex(i) = (n_real + 1i * n_imag);
    end
end


