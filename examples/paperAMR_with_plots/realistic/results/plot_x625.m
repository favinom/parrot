close all;
clear all;
clc;
figure('Position', [0, 0, 666, 500]);

hold on;

% Box-DFM
box_dfm = csvread('boxdfm/boxdfm_real_x625.csv', 1, 0);
plot(box_dfm(:, 15), box_dfm(:, 4), 'Color', 'g', 'LineWidth', 1, 'Linestyle', '-');

%% CC-DFM
box_dfm = csvread('ccdfm/ccdfm_real_x625.csv', 1, 0);
plot(box_dfm(:, 5), box_dfm(:, 1), 'Color', 'c', 'LineWidth', 1, 'Linestyle', '-');

%% CC-DFM
mpfa = csvread('ccdfm/mpfa_real_x625.csv', 1, 0);
plot(mpfa(:, 5), mpfa(:, 1), 'Color', [0 .5 0], 'LineWidth', 1, 'Linestyle', '-');


%% EDFM
edfm = csvread('edfm/edfm_real_x625.csv', 1, 0);
plot(edfm(:, 5), edfm(:, 1), 'Linestyle', '-', 'Color', 'r', 'LineWidth', 1);

% Mortar-DFM
mortar = csvread('mortardfm/mortardfm_real_x625.csv', 1, 0);
plot(mortar(:, 5), mortar(:, 1), 'Linestyle', '-', 'Color', [1 0.7 0], 'LineWidth', 1);

% format the plot
xlabel('arc length [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel('pressure [Pa]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'FontSize', 12);
legend('BoxM', 'TPFA', 'MPFA', 'EDFM', 'Flux-Mortar', 'Location', 'SouthWest');

set(gcf, 'PaperPositionMode', 'auto')
fig_pos = get(gcf, 'PaperPosition');
set(gcf, 'PaperSize', [fig_pos(3) fig_pos(4)+1])
saveas(gcf, 'real_x625', 'pdf')
