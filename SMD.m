clc
clear

% Read data from the spray_ecn.out file, skipping the 3 header lines
data = readmatrix('spray_base.out', 'FileType', 'text', 'NumHeaderLines', 4, 'Whitespace', ' \t');

% Extract columns
time_s = data(:,1); % Time in seconds (column 1)
smd_m = data(:,6);  % SMD in meters (column 6)

% Convert units
time_ms = time_s * 1000; % Convert seconds to milliseconds
smd_um = smd_m * 1e6;    % Convert meters to micrometers

%%
% Read data from the spray_ecn.out file, skipping the 3 header lines
data = readmatrix('spray_w.out', 'FileType', 'text', 'NumHeaderLines', 4, 'Whitespace', ' \t');

% Extract columns
time_s1 = data(:,1); % Time in seconds (column 1)
smd_m1 = data(:,6);  % SMD in meters (column 6)

% Convert units
time_ms1 = time_s1 * 1000; % Convert seconds to milliseconds
smd_um1 = smd_m1 * 1e6;    % Convert meters to micrometers
%%
% Read data from the spray_ecn.out file, skipping the 3 header lines
data = readmatrix('spray_s.out', 'FileType', 'text', 'NumHeaderLines', 4, 'Whitespace', ' \t');

% Extract columns
time_s2 = data(:,1); % Time in seconds (column 1)
smd_m2 = data(:,6);  % SMD in meters (column 6)

% Convert units
time_ms2 = time_s2 * 1000; % Convert seconds to milliseconds
smd_um2 = smd_m2 * 1e6;    % Convert meters to micrometers
%%
% Create the plot
figure;
plot(time_ms, smd_um, 'b-', 'LineWidth', 2, 'DisplayName', 'Base');
hold on;
plot(time_ms1, smd_um1, 'r-', 'LineWidth', 2, 'DisplayName', 'Thermal-Breakup-W');
hold on;
plot(time_ms2, smd_um2, 'g-', 'LineWidth', 2, 'DisplayName', 'Thermal-Breakup-S');
hold on;

% Add labels and title
xlabel('Time (ms)');
ylabel('SMD (\mum)'); % \mum for micrometers symbol
title('Sauter Mean Diameter vs. Time');
grid on;
legend('show');

% Set font size for better readability
set(gca, 'FontSize', 12);

