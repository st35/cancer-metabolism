clc;
clear all;

P1 = logspace(0, 3, 50);
P2 = logspace(-2, 3, 50);

N_Metabolites = 28;
N_Fluxes = 32;

F0 = zeros(length(P1), length(P2));
F1 = zeros(length(P1), length(P2));

Input = zeros(2, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.160 - x0(24); % M_ADP
x0(26) = 1.0; % C_ATP
x0(28) = 0.001; % C_NADH

for i = length(P1):-1:1
    for j = 1:length(P2)
        [num2str(i) ' ' num2str(j)]
        Input(1) = P1(i);
        Input(2) = P2(j);
        [t, x] = ode23tb(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 5000000.0], x0);
        x0 = x(end, :)';
        F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

        if sum(isnan(x0)) > 0
            'Dying is easy. Integrating this is hard. - James Wilson.'
        end

        F0(i, j) = F(30); % NAD use
        F1(i, j) = F(31); % NAD:NADH ratio
    end
end

figure;
imagesc(log10(P2), log10(P1), F0);
colormap redbluecmap;
set(gca, 'FontSize', 16, 'Box', 'on', 'YDir', 'normal');
caxis([150 500]);
colorbar;
xlabel('Factor change in PDH activity');
ylabel('Factor change in LDH activity');
title('NAD^{+} use flux', 'FontWeight', 'normal');
xticks(-2:3);
set(gca, 'xticklabel', {'10^{-2}', '10^{-1}', '10^0', '10^1', '10^2', '10^3'});
yticks(0:3);
set(gca, 'yticklabel', {'10^0', '10^1', '10^2', '10^3'});

figure;
imagesc(log10(P2), log10(P1), log2(F1));
colormap redbluecmap;
set(gca, 'FontSize', 16, 'Box', 'on', 'YDir', 'normal');
caxis([4 9]);
colorbar;
xlabel('Factor change in PDH activity');
ylabel('Factor change in LDH activity');
title('NAD^{+}:NADH ratio (log2)', 'FontWeight', 'normal');
xticks(-2:3);
set(gca, 'xticklabel', {'10^{-2}', '10^{-1}', '10^0', '10^1', '10^2', '10^3'});
yticks(0:3);
set(gca, 'yticklabel', {'10^0', '10^1', '10^2', '10^3'});