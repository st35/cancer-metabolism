clc;
clear all;

P1 = logspace(-3, 3, 100);
P2 = logspace(-3, 3, 100);

N_Metabolites = 31;
N_Fluxes = 38;

F0 = zeros(length(P1), length(P2));
F1 = zeros(length(P1), length(P2));

Input = zeros(2, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.16*4.0 - x0(24); % M_ADP
x0(26) = 1.0; % C_ATP

for i = 1:length(P1)
    x0 = zeros(N_Metabolites, 1);
    x0(21) = 4.0e-2*1e2; % M_OAA
    x0(24) = 900.0e-4; % M_ATP
    x0(25) = 4.16*4.0 - x0(24); % M_ADP
    x0(26) = 1.0; % C_ATP
    for j = 1:length(P2)
        % x0 = zeros(N_Metabolites, 1);
        % x0(21) = 4.0e-2*1e2; % M_OAA
        % x0(24) = 900.0e-4; % M_ATP
        % x0(25) = 225.0e-4; % M_ADP
        % x0(26) = 1.0; % C_ATP
        [num2str(i) ' ' num2str(j)]
        Input(1) = P1(i);
        Input(2) = P2(j);
        [t, x] = ode23tb(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 50000.0], x0);
        x0 = x(end, :)';
        F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

        if sum(isnan(x0)) > 0
            'Dying is easy. Integrating this is hard. - James Wilson.'
        end

        F0(i, j) = F(20);
        F1(i, j) = F(38);
    end
end

figure;
% hold on;
imagesc(log10(P2), log10(P1), F0);
set(gca, 'FontSize', 16, 'Box', 'on', 'YDir', 'normal');
colorbar;
colormap redbluecmap;
caxis([min(min(F0)) max(max(F0))]);
title('TCA cycle flux', 'FontWeight', 'normal');
xlabel('Factor change in ACC activity');
ylabel('Glutamine conc.');
xticks(-3:3);
yticks(-3:3);
set(gca, 'xticklabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}', '10^{3}',})
set(gca, 'yticklabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}', '10^{3}',})

figure;
% hold on;
imagesc(log10(P2), log10(P1), F1);
set(gca, 'FontSize', 16, 'Box', 'on', 'YDir', 'normal');
colorbar;
colormap redbluecmap;
caxis([min(min(F1)) max(max(F1))]);
title('Fatty acids syn. flux', 'FontWeight', 'normal');
xlabel('Factor change in ACC activity');
ylabel('Glutamine conc.');
xticks(-3:3);
yticks(-3:3);
set(gca, 'xticklabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}', '10^{3}',})
set(gca, 'yticklabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}', '10^{3}',})