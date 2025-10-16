% Satellite Constellation Simulation and Graphics - Matlab CODE
clear; clc; close all;

% Define the earth and sim parameters
mu = 3.986e5;                  % Earth's gravitational parameter [km^3/s^2]

Radius = 6371;                                          % Earth Radius [km]

altitude = input("Enter the esimated altitude of the constellation [km]: ");                                      
                                                     % above the earth [km]

a = Radius + altitude;                               % Semi-Major Axis [km]

e = 0.001;                                            % Near circular orbit

i = deg2rad(55);                                        % inclination [rad]

omega = deg2rad(0);                            % Argument of perigree [rad]


% Constellation Parameters
satnum = input('Enter the number of satellites in the constellation: ');

% Distribute the satellites in RAAN (and optionally phase) uniformally:
RAANs = linspace(0,2*pi,satnum+1);
RAANs(end) = []; % this is to remove the duplicate 360 deg 
%Offset initial mean anomaly for diversity
phase_offsets = linspace(0, 2*pi, satnum+1);
phase_offsets(end) = [];

% Orbital Period
T = 2*pi*sqrt(a^3/mu);
tspan = linspace(0, T, 1000); % the period divided into steps

% Figure setup
figure;
hold on;
axis equal;
grid on;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Satellite constellation Simulation');

% Plot the earth
[xe,ye,ze] = sphere(50);
surf(xe*Radius, ye*Radius, ze*Radius, 'FaceColor','blue', 'EdgeColor','none','FaceAlpha',0.3);

% Loop over each satellite in the constellation
for s = 1:satnum
    RAAN = RAANs(s);
    M0 = phase_offsets(s); % Initial mean anomaly offset

    % Initialize the matrix to store positions
    positions = zeros(3, length(tspan));

    % Compute the satellites position at each time step
    for k = 1:length(tspan)
        
        % Mean anomaly involves lineaerly with time. 
        M = M0 + 2*pi*tspan(k)/T;

        % Solve for the eccentric anamoly, E, using Kepler's equation with iteration
        E = M; % Initial guess
        for iter = 1:100
            E = E-(E-e*sin(E)-M)/(1-e*cos(E));
        end
        
        % Compute the true anomaly, nu
        nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1+e)*cos(E/2));

        % Perifocal coordinates (Orbital plane)
        perifocal = (a*(1-e^2))/(1+e*cos(nu)) * [cos(nu); sin(nu); 0];

        % Compute the rotation matrix from the perifocal to the ECI coordinates
        R3_W = [cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];

        R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    
        R3_w = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];

        Q = (R3_W')*(R1_i')*(R3_w');

        % Convert to ECI coordinates
        r_eci = Q*perifocal;
        positions(:,k) = r_eci;
    end
    
    %plot the orbit trajectory and mark the starting position
    plot3(positions(1,:), positions(2,:), positions(3,:), 'LineWidth', 1.5);
    plot3(positions(1,1), positions(2,1), positions(3,1),'o', 'MarkerSize',6, 'MarkerFaceColor','r');
end

view(3); % 3D view
hold off;


