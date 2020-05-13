close all;
clear all;
clc;
figure('Position', [0, 0, 666, 500]);

ishoriz=2;
%hold on;

filenames{1,1}='boxdfm/boxdfm_real_y500.csv';
filenames{2,1}='ccdfm/mpfa_real_y500.csv';
filenames{3,1}='ccdfm/ccdfm_real_y500.csv';
filenames{4,1}='edfm/edfm_real_y500.csv';
filenames{5,1}='mortardfm/mortardfm_real_y500.csv';

filenames{1,2}='boxdfm/boxdfm_real_x625.csv';
filenames{2,2}='ccdfm/mpfa_real_x625.csv';
filenames{3,2}='ccdfm/ccdfm_real_x625.csv';
filenames{4,2}='edfm/edfm_real_x625.csv';
filenames{5,2}='mortardfm/mortardfm_real_x625.csv';

myfilename{1}='../horizontal_line_1_12_0.csv';
myfilename{2}='../vertical_line_1_12_0.csv';

% Box-DFM
for i=1:5
    res{i}=csvread(filenames{i,ishoriz}, 1, 0);
end

mylegend{1}='Box';
mylegend{2}='TPFA';
mylegend{3}='MPFA';
mylegend{4}='EDFM';
mylegend{5}='Flux-Mortar';

mycolor{1}='g';
mycolor{2}='c';
mycolor{3}=[0 .5 0];
mycolor{4}='r';
mycolor{5}=[1 0.7 0];

myres = csvread(myfilename{ishoriz}, 1);


for i=1:5
subplot(2,3,i)
if (i==1)
    plot(res{i}(:, 15), res{i}(:, 4), 'Color', mycolor{i}, 'LineWidth', 1, 'Linestyle', '-');
else
    if ishoriz==1
        plot(res{i}(:, 4), res{i}(:, 1), 'Color', mycolor{i}, 'LineWidth', 1, 'Linestyle', '-');
    else
        plot(res{i}(:, 5), res{i}(:, 1), 'Color', mycolor{i}, 'LineWidth', 1, 'Linestyle', '-');
    end
end    
hold on
if (ishoriz==1)
    plot(myres(:,4),myres(:,1), 'Color', 'k', 'LineWidth', 1, 'Linestyle', '-');
else
    plot(myres(:,5),myres(:,1), 'Color', 'k', 'LineWidth', 1, 'Linestyle', '-');
end
title(mylegend{i},'FontSize', 10,'FontWeight', 'normal')
xlabel('arc length [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel('pressure [Pa]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'FontSize', 12);

end
%subplot(2,3,6)

set(gcf, 'PaperPositionMode', 'auto')
fig_pos = get(gcf, 'PaperPosition');
set(gcf, 'PaperSize', [fig_pos(3) fig_pos(4)+1])
saveas(gcf, 'real_y500', 'pdf')


return

% format the plot


