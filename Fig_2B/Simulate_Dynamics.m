clc;
clear all;

P1 = logspace(-4, 0, 100);
P2 = logspace(-4, -1, 100);

N_Metabolites = 26;
N_Fluxes = 31;

F0 = zeros(length(P1), length(P2));

Input = zeros(2, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.16*4.0 - x0(24); % M_ADP
x0(26) = 1.0; % C_ATP

for i = length(P1):-1:1
    x0 = zeros(N_Metabolites, 1);
    x0(21) = 4.0e-2*1e2; % M_OAA
    x0(24) = 900.0e-4; % M_ATP
    x0(25) = 4.16*4.0 - x0(24); % M_ADP
    x0(26) = 1.0; % C_ATP
    for j = length(P2):-1:1
        [num2str(i) ' ' num2str(j)]
        Input(1) = P1(i);
        Input(2) = P2(j);
        [t, x] = ode23tb(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 50000.0], x0);
        x0 = x(end, :)';
        F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

        if sum(isnan(x0)) > 0
            'Dying is easy. Integrating this is hard. - James Wilson.'
        end

        F0(i, j) = F(14);
    end
end

figure;
imagesc(log10(P2), log10(P1), F0);
colormap redbluecmap;
set(gca, 'FontSize', 16, 'Box', 'on', 'YDir', 'normal');
colorbar;
caxis([min(min(F0)) max(max(F0))]);
title('Lactate production flux', 'FontWeight', 'normal');
xlabel('Factor change in rate of ATP use');
ylabel('O_{2} concentration');
set(gca, 'xticklabel', {'10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}'});
set(gca, 'yticklabel', {'10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{0}'});