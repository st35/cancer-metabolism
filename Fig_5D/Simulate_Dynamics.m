clc;
clear all;

P1 = logspace(-3, 1, 100);
P2 = logspace(-2, 0, 100);

N_Metabolites = 31;
N_Fluxes = 39;

F0 = zeros(length(P1), length(P2));

Input = zeros(2, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(22) = 1.0; % M_NADH
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.0*4.16 - x0(25); % M_ADP
x0(26) = 1.0; % C_ATP

for i = length(P1):-1:1
    x0 = zeros(N_Metabolites, 1);
    x0(21) = 4.0e-2*1e2; % M_OAA
    % x0(22) = 1.0; % M_NADH
    x0(24) = 900.0e-4; % M_ATP
    x0(25) = 4.0*4.16 - x0(25); % M_ADP
    x0(26) = 1.0; % C_ATP
    for j = 1:length(P2)
        [num2str(i) ' ' num2str(j)]
        Input(1) = P1(i);
        Input(2) = P2(j);
        [t, x] = ode23s(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 50000.0], x0);
        x0 = x(end, :)';
        F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

        if sum(isnan(x0)) > 0
            'Dying is easy. Integrating this is hard. - James Wilson.'
        end

        F0(i, j) = F(28);
    end
end

figure;
imagesc(log10(P2), log10(P1), F0);
colormap redbluecmap;
set(gca, 'FontSize', 16, 'Box', 'on', 'YDir', 'normal');
colorbar;
caxis([min(min(F0)) max(max(F0))]);
xlabel('O_{2} conc.');
ylabel({'Factor change in ME2 activity'});
title('ATP use flux', 'FontWeight', 'normal');
xticks(-2:0);
yticks(-3:1);
set(gca, 'xticklabel', {'10^{-2}', '10^{-1}', '10^{0}'});
set(gca, 'yticklabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^{0}', '10^{1}'});