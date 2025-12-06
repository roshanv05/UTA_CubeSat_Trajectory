%Orbit result plot

% reads and plots orbit from orbit_results.csv
% Date: 12/05/25
% Author: Karthik
% Affiliation: UTA AIAA
% make sure the output from the 2body_orbital_dynamics.cpp is in the same file directory as the activated matlab.
% ============================================================

clear; clc; close all;


%% === READ CSV output from C++ ===
filename = 'orbit_results';
% M will be 1 x 20 (one data row, 20 columns), skipping the header
M = readmatrix(filename, 'NumHeaderLines', 1);

% Column mapping (must match the C++ header order):
%  1: a_leo_km
%  2: e_leo
%  3: r_peri_leo_km
%  4: r_apo_leo_km
%  5: a_sso_km
%  6: e_sso
%  7: inc_sso_deg
%  8: r_peri_sso_km
%  9: r_apo_sso_km
% 10: a_trans_km
% 11: e_trans
% 12: r_peri_trans_km
% 13: r_apo_trans_km
% 14: dv1_kms
% 15: dv2_kms
% 16: dv_total_kms
% 17: m0_kg
% 18: mf_kg
% 19: mprop_kg
% 20: Itot_Ns

a_leo      = M(1,1);     % [km]
e_leo      = M(1,2);
a_sso      = M(1,5);     % [km]
e_sso      = M(1,6);
inc_sso    = deg2rad(M(1,7)); % [rad]
a_trans    = M(1,10);    % [km]
e_trans    = M(1,11);

% You can change these if you want multiple sats, different RAAN, etc.
RAAN_leo   = 0;                % [rad]
RAAN_sso   = 0;                % [rad]
RAAN_trans = 0;                % [rad]
omega_leo   = 0;               % argument of perigee [rad]
omega_sso   = 0;               % [rad]
omega_trans = 0;               % [rad]

% Inclinations:
i_leo   = 0;          % equatorial LEO (2D-ish)
i_sso   = inc_sso;    % from C++ sun-synchronous calc
i_trans = i_leo;      % same plane as LEO for this simple model


%% === CONSTANTS ===
mu     = 3.986e5;    % Earth's gravitational parameter [km^3/s^2]
Radius = 6371;       % Earth radius [km]

% Time discretization: one period for each orbit
T_leo   = 2*pi*sqrt(a_leo^3 / mu);
T_sso   = 2*pi*sqrt(a_sso^3 / mu);
T_trans = 2*pi*sqrt(a_trans^3 / mu);

Npts = 1000;   % number of points per orbit
t_leo   = linspace(0, T_leo,   Npts);
t_sso   = linspace(0, T_sso,   Npts);
t_trans = linspace(0, T_trans, Npts);


%% === HELPER FUNCTION FOR KEPLER + ORBIT PROPAGATION ===
propagate_orbit = @(a, e, inc, RAAN, omega, tspan) ...
    local_propagate_orbit(a, e, inc, RAAN, omega, tspan, mu);


%% === PROPAGATE ORBITS ===
% LEO
positions_leo = propagate_orbit(a_leo, e_leo, i_leo, RAAN_leo, omega_leo, t_leo);

% SSO
positions_sso = propagate_orbit(a_sso, e_sso, i_sso, RAAN_sso, omega_sso, t_sso);

% Transfer orbit
positions_trans = propagate_orbit(a_trans, e_trans, i_trans, RAAN_trans, omega_trans, t_trans);


%% === PLOTTING (similar style to your constellation code) ===
figure; hold on; axis equal; grid on;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('LEO, SSO, and Transfer Orbits');

% Plot Earth
[xe, ye, ze] = sphere(50);
surf(xe*Radius, ye*Radius, ze*Radius, ...
    'FaceColor','blue', 'EdgeColor','none','FaceAlpha',0.3);

% Plot LEO orbit
plot3(positions_leo(1,:), positions_leo(2,:), positions_leo(3,:), ...
      'LineWidth', 1.5);
plot3(positions_leo(1,1), positions_leo(2,1), positions_leo(3,1), ...
      'o', 'MarkerSize', 6, 'MarkerFaceColor','r');

% Plot SSO orbit
plot3(positions_sso(1,:), positions_sso(2,:), positions_sso(3,:), ...
      'LineWidth', 1.5);
plot3(positions_sso(1,1), positions_sso(2,1), positions_sso(3,1), ...
      'o', 'MarkerSize', 6, 'MarkerFaceColor','g');

% Plot transfer orbit
plot3(positions_trans(1,:), positions_trans(2,:), positions_trans(3,:), ...
      'LineWidth', 1.5, 'LineStyle','--');
plot3(positions_trans(1,1), positions_trans(2,1), positions_trans(3,1), ...
      'o', 'MarkerSize', 6, 'MarkerFaceColor','k');

legend('Earth', 'LEO orbit', 'LEO start', ...
       'SSO orbit', 'SSO start', ...
       'Transfer orbit', 'Transfer start');

view(3);
hold off;



%% === LOCAL FUNCTION: PROPAGATE ORBIT USING KEPLER'S EQUATION ===
function positions = local_propagate_orbit(a, e, inc, RAAN, omega, tspan, mu)
    % Returns positions (3 x N) in ECI for one orbit using Kepler's equation

    % Orbital period
    T = 2*pi*sqrt(a^3/mu);

    % Preallocate
    N = length(tspan);
    positions = zeros(3, N);

    % Loop over time
    for k = 1:N
        t = tspan(k);

        % Mean anomaly evolves linearly with time
        M = 2*pi*t / T;

        % Solve Kepler's equation for eccentric anomaly E using iteration
        E = M;  % initial guess
        for iter = 1:50
            E = E - (E - e*sin(E) - M) / (1 - e*cos(E));
        end

        % True anomaly
        nu = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );

        % Perifocal radius
        r = (a*(1 - e^2)) / (1 + e*cos(nu));

        % Perifocal coordinates
        perifocal = r * [cos(nu); sin(nu); 0];

        % Rotation matrices (same style as your constellation code)
        R3_W = [ cos(RAAN)  sin(RAAN)  0;
                -sin(RAAN)  cos(RAAN)  0;
                 0          0          1];

        R1_i = [ 1      0           0;
                 0  cos(inc)   sin(inc);
                 0 -sin(inc)   cos(inc)];

        R3_w = [ cos(omega)  sin(omega)  0;
                -sin(omega)  cos(omega)  0;
                 0           0           1];

        Q = (R3_W') * (R1_i') * (R3_w');

        % Convert to ECI coordinates
        r_eci = Q * perifocal;

        positions(:,k) = r_eci;
    end
end
