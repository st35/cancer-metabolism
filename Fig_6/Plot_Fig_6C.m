clc;
clear all;

D = importdata('data/GSE75299_PatientFPKM.txt');
data = D.data;
L = D.textdata;
Sample_Names = strsplit(L{1});
Gene_Names = {L{2:end}};

PC_index = find(strcmp(Gene_Names, 'PC'));
PC = log2(data(PC_index, :));
PC = (PC - mean(PC)) / std(PC);

MYC_index = find(strcmp(Gene_Names, 'MYC'));
MYC = log2(data(MYC_index, :));
MYC = (MYC - mean(MYC)) / std(MYC);

map = brewermap(3, 'Set1');
figure;

Pat_Index = find(contains(Sample_Names, 'Pt1'));

subplot(2, 3, 1);
hold on;
for i = 1:length(Pat_Index)
    if i == 1
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(2, :), 'filled');
    elseif i == 2
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), '*', 'linewidth', 1.5);
    end
end
title('Patient 1', 'FontWeight', 'normal');
set(gca, 'FontSize', 16, 'Box', 'on');
xlabel('MYC');
ylabel('PC');

Pat_Index = find(contains(Sample_Names, 'Pt3'));

subplot(2, 3, 2);
hold on;
for i = 1:length(Pat_Index)
    if i == 1
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(2, :), 'filled');
    elseif i == 2
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), '*', 'linewidth', 1.5);
    end
end
title('Patient 3', 'FontWeight', 'normal');
set(gca, 'FontSize', 16, 'Box', 'on');
xlabel('MYC');
ylabel('PC');

Pat_Index = find(contains(Sample_Names, 'Pt4'));

subplot(2, 3, 3);
hold on;
for i = 1:length(Pat_Index)
    if i == 1
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(2, :), 'filled');
    elseif i == 2
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), '*', 'linewidth', 1.5);
    elseif i == 3
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), 'x', 'linewidth', 1.5);
    end
end
title('Patient 4', 'FontWeight', 'normal');
set(gca, 'FontSize', 16, 'Box', 'on');
xlabel('MYC');
ylabel('PC');

Pat_Index = find(contains(Sample_Names, 'Pt6'));

subplot(2, 3, 4);
hold on;
for i = 1:length(Pat_Index)
    if i == 1
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(2, :), 'filled');
    elseif i == 2
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), '*', 'linewidth', 1.5);
    elseif i == 3
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), 'x', 'linewidth', 1.5);
    else
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), '+', 'linewidth', 1.5);
    end
end
title('Patient 6', 'FontWeight', 'normal');
set(gca, 'FontSize', 16, 'Box', 'on');
xlabel('MYC');
ylabel('PC');

Pat_Index = find(contains(Sample_Names, 'Pt7'));

subplot(2, 3, 5);
hold on;
for i = 1:length(Pat_Index)
    if i == 1
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(2, :), 'filled');
    elseif i == 2
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), '*', 'linewidth', 1.5);
    elseif i == 3
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), 'x', 'linewidth', 1.5);
    else
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), '+', 'linewidth', 1.5);
    end
end
title('Patient 7', 'FontWeight', 'normal');
set(gca, 'FontSize', 16, 'Box', 'on');
xlabel('MYC');
ylabel('PC');

Pat_Index = find(contains(Sample_Names, 'Pt8'));

subplot(2, 3, 6);
hold on;
for i = 1:length(Pat_Index)
    if i == 1
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(2, :), 'filled');
    elseif i == 2
        scatter(MYC(Pat_Index(i)), PC(Pat_Index(i)), 160, map(1, :), '*', 'linewidth', 1.5);
    end
end
title('Patient 8', 'FontWeight', 'normal');
set(gca, 'FontSize', 16, 'Box', 'on');
xlabel('MYC');
ylabel('PC');