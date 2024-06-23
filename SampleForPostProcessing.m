clear; clc; close all; % MATLAB Environment cleanup

%PROGRAM FOR SYNTHETIZING AND PRODUCING NET PLOTS FOR ALL
%PLASMA PROPERTIES' COMPUTATIONS.

% List of file names (modify according to file names)
fileNames = {'50kW.txt', '100kW.txt', '150kW.txt', '200kW.txt', '250kW.txt', '300kW.txt'};

% Initialize cell array to store data from each file
data = cell(length(fileNames), 1);

% Load data from each file
for i = 1:length(fileNames)
    data{i} = readmatrix(fileNames{i});
end

% Extract the fixed x-axis data from the first file (assume it's the first column)
Z = data{1}(:, 1);

% Create a new figure for plotting
figure(4);
clf
set(gcf, 'Position', get(0, 'Screensize'));

%Fixed parameters
mDot = 80;
En = 0.1;
SHrat = 5/4;
alpha_exc = 1;

%Normalization parameters


%Control variable for legends
P=50;

% Plot data from each file
for i = 1:length(fileNames)
    % Extract y-axis data (assume it's the second column)
    ne = data{i}(:, 2);
    nn = data{i}(:, 3);
    Te = data{i}(:, 4);
    uzi = data{i}(:, 5);
    ASp = data{i}(:, 6);
    mFlow = data{i}(:, 7);

    % Plot of values against normalized system length
    subplot(2, 3, 1);
    hold on;
    plot(Z, ne, 'DisplayName', sprintf('%.2f kW', P));
    title('Electron Density');
    xlabel('Normalized Length');
    ylabel('Normalized Electron Density');
    legend show;
    hold off;

    %Neutral Density
    subplot(2, 3, 2);
    hold on;
    plot(Z, nn, 'DisplayName', sprintf('%.2f kW', P));
    title('Neutral Density');
    xlabel('Normalized Length');
    ylabel('Normalized Neutral Density');
    legend show;
    hold off;

    %Electron Temperature
    subplot(2, 3, 3);
    hold on;
    plot(Z, Te, 'DisplayName', sprintf('%.2f kW', P));
    title('Electron Temperature');
    xlabel('Normalized Length');
    ylabel('Normalized Electron Temperature');
    legend show;
    hold off;

    %Ion Velocity
    subplot(2, 3, 4);
    hold on;
    plot(Z, uzi, 'DisplayName', sprintf('%.2f kW', P));
    title('Ion Velocities');
    xlabel('Normalized Length');
    ylabel('Normalized Ion Velocity');
    legend show;
    hold off;

    %Ionization Source Amplitude
    subplot(2, 3, 5);
    hold on;
    plot(Z, ASp, 'DisplayName', sprintf('%.2f kW', P));
    title('Plasma Ionization');
    xlabel('Normalized Length');
    ylabel('Normalized Plasma Ionization');
    legend show;
    hold off;

    %Mass Flow
    subplot(2, 3, 6);
    hold on;
    plot(Z, mFlow, 'DisplayName', sprintf('%.2f kW', P));
    title('Mass Flow');
    xlabel('Normalized Length');
    ylabel('Normalized Mass Flow');
    legend show;
    hold off;

    P=P+50; %Variables were computed in steps of 50kW.
end

sgtitle(sprintf("Plasma Dynamics of a DFD for Low-Power Inputs\n mDot = %0.f mg/s | En = %.2f eV | SHrat = %.2f | alpha = %1.f", mDot, En, SHrat, alpha_exc));