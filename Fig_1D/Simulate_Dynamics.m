clc;
clear all;

P1 = logspace(-3, 0, 100);
P2 = logspace(-2, 2, 100);

N_Metabolites = 21;
N_Fluxes = 24;

F0 = zeros(length(P1), length(P2));

Input = zeros(2, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA

for i = length(P1):-1:1
    i
    for j = 1:length(P2)
        Input(1) = P1(i);
        Input(2) = P2(j);
        [t, x] = ode23s(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 5000000.0], x0);
        x0 = x(end, :)';
        F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

        if sum(isnan(x0)) > 0
            'Dying is easy. Integrating this is hard. - James Wilson.'
        end

        F0(i, j) = F(8);
    end
end

figure;
imagesc(log10(P2), log10(P1), F0);
colormap redbluecmap;
set(gca, 'FontSize', 16, 'Box', 'on', 'YDir', 'normal');
colorbar;
caxis([min(min(F0)) max(max(F0))]);
title('Phospholipids syn. flux', 'FontWeight', 'normal');
xlabel('Factor change in PDH activity');
ylabel('Factor change in LDH activity');
set(gca, 'xticklabel', {'10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}'});
yticks(-3:0);
set(gca, 'yticklabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^{0}'});