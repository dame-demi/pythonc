clc
clear

% Read data from the spray_ecn.out file, skipping the 3 header lines
data = readmatrix('spray_ecn_nothermal.out', 'FileType', 'text', 'NumHeaderLines', 3, 'Whitespace', ' \t');
L = [   0.0000000e+00   0.0000000e+00   0.0000000e+00
   5.0000000e+01   4.0280000e+00   1.7389698e+00
   1.0000000e+02   9.8760000e+00   1.5040691e+00
   1.5000000e+02   1.4752000e+01   1.2740863e+00
   2.0000000e+02   1.8996000e+01   1.1443706e+00
   2.5000000e+02   2.3132000e+01   1.0564923e+00
   3.0000000e+02   2.6912000e+01   1.0057117e+00
   3.5000000e+02   3.0264000e+01   9.5869912e-01
   4.0000000e+02   3.3288000e+01   9.7992653e-01
   4.5000000e+02   3.6184000e+01   9.6733862e-01
   5.0000000e+02   3.8708000e+01   9.1582531e-01
   5.5000000e+02   4.0948000e+01   9.6233882e-01
   6.0000000e+02   4.2932000e+01   1.0342998e+00
   6.5000000e+02   4.4896000e+01   9.4593023e-01
   7.0000000e+02   4.6684000e+01   8.9718671e-01
   7.5000000e+02   4.8412000e+01   8.9836296e-01
   8.0000000e+02   5.0096000e+01   9.6849574e-01
   8.5000000e+02   5.1716000e+01   1.0763568e+00
   9.0000000e+02   5.3272000e+01   1.2754670e+00
   9.5000000e+02   5.4780000e+01   1.3753545e+00
   1.0000000e+03   5.6212000e+01   1.4313127e+00
   1.0500000e+03   5.7628000e+01   1.5200053e+00
   1.1000000e+03   5.9016000e+01   1.6242364e+00
   1.1500000e+03   6.0316000e+01   1.7376260e+00
   1.2000000e+03   6.1544000e+01   1.8829934e+00
   1.2500000e+03   6.2696000e+01   2.0132521e+00
   1.3000000e+03   6.3744000e+01   2.1234086e+00
   1.3500000e+03   6.4752000e+01   2.1464613e+00
   1.4000000e+03   6.5840000e+01   2.1402804e+00
   1.4500000e+03   6.6772000e+01   2.1967285e+00
   1.5000000e+03   6.7624000e+01   2.2089418e+00
   1.5500000e+03   6.8400000e+01   2.2676860e+00
   1.6000000e+03   6.9172000e+01   2.2467790e+00
   1.6500000e+03   6.9804000e+01   2.2561879e+00
   1.7000000e+03   7.0384000e+01   2.1754411e+00
   1.7500000e+03   7.0916000e+01   2.0611026e+00
   1.8000000e+03   7.1368000e+01   1.9828706e+00
   1.8500000e+03   7.1712000e+01   1.8906761e+00
   1.9000000e+03   7.1980000e+01   1.8063222e+00
   1.9500000e+03   7.2284000e+01   1.6647354e+00
   2.0000000e+03   7.1140000e+01   1.0273558e+01];
% Extract time (column 1)
time = data(:, 1);
L(:,1) = L(:,1)/1000;

% Define column indices for Spray_Penet97 (liquid) and Vapor_Penet for each nozzle
liquid_cols = [6, 12, 18, 24, 30, 36, 42, 48]; % Spray_Penet97 columns
vapor_cols = [8, 14, 20, 26, 32, 38, 44, 50];   % Vapor_Penet columns

% Extract liquid and vapor penetration data
liquid_penet = data(:, liquid_cols);
vapor_penet = data(:, vapor_cols);

% Compute average penetration across the eight nozzles
avg_liquid_penet = mean(liquid_penet, 2);
avg_vapor_penet = mean(vapor_penet, 2);

% Convert time from seconds to milliseconds for better readability
time_ms = time * 1e3;

% Convert penetration from meters to millimeters
avg_liquid_penet_mm = avg_liquid_penet * 1e3;
avg_vapor_penet_mm = avg_vapor_penet * 1e3;

%%
% Read data from the spray_ecn.out file, skipping the 3 header lines
data1 = readmatrix('spray_ecn_thermal1.out', 'FileType', 'text', 'NumHeaderLines', 3, 'Whitespace', ' \t');

% Extract time (column 1)
time = data1(:, 1);

% Define column indices for Spray_Penet97 (liquid) and Vapor_Penet for each nozzle
liquid_cols = [6, 12, 18, 24, 30, 36, 42, 48]; % Spray_Penet97 columns
vapor_cols = [8, 14, 20, 26, 32, 38, 44, 50];   % Vapor_Penet columns

% Extract liquid and vapor penetration data
liquid_penet1 = data1(:, liquid_cols);
vapor_penet1 = data1(:, vapor_cols);

% Compute average penetration across the eight nozzles
avg_liquid_penet1 = mean(liquid_penet1, 2);
avg_vapor_penet1 = mean(vapor_penet1, 2);

% Convert time from seconds to milliseconds for better readability
time_ms1 = time * 1e3;

% Convert penetration from meters to millimeters
avg_liquid_penet_mm1 = avg_liquid_penet1 * 1e3;
avg_vapor_penet_mm1 = avg_vapor_penet1 * 1e3;
%%
% Read data from the spray_ecn.out file, skipping the 3 header lines
data2 = readmatrix('spray_ecn_thermal2.out', 'FileType', 'text', 'NumHeaderLines', 3, 'Whitespace', ' \t');

% Extract time (column 1)
time = data2(:, 1);

% Define column indices for Spray_Penet97 (liquid) and Vapor_Penet for each nozzle
liquid_cols = [6, 12, 18, 24, 30, 36, 42, 48]; % Spray_Penet97 columns
vapor_cols = [8, 14, 20, 26, 32, 38, 44, 50];   % Vapor_Penet columns

% Extract liquid and vapor penetration data
liquid_penet2 = data2(:, liquid_cols);
vapor_penet2 = data2(:, vapor_cols);

% Compute average penetration across the eight nozzles
avg_liquid_penet2 = mean(liquid_penet2, 2);
avg_vapor_penet2 = mean(vapor_penet2, 2);

% Convert time from seconds to milliseconds for better readability
time_ms2 = time * 1e3;

% Convert penetration from meters to millimeters
avg_liquid_penet_mm2 = avg_liquid_penet2 * 1e3;
avg_vapor_penet_mm2 = avg_vapor_penet2 * 1e3;

%%
% Create the plot
figure(1);
plot(time_ms, avg_liquid_penet_mm, 'b-', 'LineWidth', 2, 'DisplayName', 'NoThermal');
hold on;
plot(time_ms1, avg_liquid_penet_mm1, 'r-', 'LineWidth', 2, 'DisplayName', 'Thermal-Breakup-W');
hold on;
plot(time_ms2, avg_liquid_penet_mm2, 'g-', 'LineWidth', 2, 'DisplayName', 'Thermal-Breakup-S');
hold on;
errorbar(Z(:,1), Z(:,2), Z(:,3), 'k-', 'LineWidth', .2, 'DisplayName', 'Experimental with Error');
plot(Z(:,1), Z(:,2), 'k.', 'MarkerSize', 10, 'DisplayName', 'Experimental Points');
% Add labels and title
xlabel('Time (ms)');
ylabel('Penetration (mm)');
title('Average Liquid Penetration vs. Time');
grid on;
legend('show');

% Set font size for better readability
set(gca, 'FontSize', 12);

%%
% Create the plot
figure(2);
plot(time_ms, avg_vapor_penet_mm, 'b-', 'LineWidth', 2, 'DisplayName', 'NoThermal');
hold on;
plot(time_ms1, avg_vapor_penet_mm1, 'r-', 'LineWidth', 2, 'DisplayName', 'Thermal-Breakup-W');
hold on;
plot(time_ms2, avg_vapor_penet_mm2, 'g-', 'LineWidth', 2, 'DisplayName', 'Thermal-Breakup-S');

% Add labels and title
xlabel('Time (ms)');
ylabel('Penetration (mm)');
title('Average Vapor Penetration vs. Time');
grid on;
legend('show');

% Set font size for better readability
set(gca, 'FontSize', 12);