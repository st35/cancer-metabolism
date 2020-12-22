clc;
clear all;

D = importdata('data/Jia_Paudel_Front_Onc_Subclones_PC.log');
D = D.data;
D = log2(D);
D = (D - mean(D)) / std(D);

SC01 = reshape(D(1:9), [3, 3]);
SC07 = reshape(D(10:18), [3, 3]);
SC10 = reshape(D(19:27), [3, 3]);

Days = [0 3 8];

map = brewermap(3, 'Dark2');

figure;
hold on;
h1 = plot(Days, mean(SC01), '-o', 'color', map(1, :), 'linewidth', 2);
errorbar(Days, mean(SC01), std(SC01), '.', 'color', map(1, :), 'linewidth', 2);
h2 = plot(Days, mean(SC07), '-o', 'color', map(2, :), 'linewidth', 2);
errorbar(Days, mean(SC07), std(SC07), '.', 'color', map(2, :), 'linewidth', 2);
h3 = plot(Days, mean(SC10), '-o', 'color', map(3, :), 'linewidth', 2);
errorbar(Days, mean(SC10), std(SC10), '.', 'color', map(3, :), 'linewidth', 2);
set(gca, 'FontSize', 16, 'Box', 'on');
legend([h1 h2 h3], {'SC01', 'SC07', 'SC10'}, 'location', 'southwest');
xlabel('Time (Days)');
ylabel('PC');