clc;
clear all;

P1 = logspace(-3, 1, 500);

N_Metabolites = 26;
N_Fluxes = 30;

F0 = zeros(length(P1), 1);
F1 = zeros(length(P1), 1);

Input = zeros(1, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.16*4.0 - x0(24); % M_ADP
x0(26) = 1.0; % C_ATP

for i = 1:length(P1)
    i
    Input(1) = P1(i);
    [t, x] = ode23tb(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 50000.0], x0);
    x0 = x(end, :)';
    F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

    if sum(isnan(x0)) > 0
        'Dying is easy. Integrating this is hard. - James Wilson.'
    end

    F0(i) = F(8);
    F1(i) = F(28);
end

map = brewermap(3, 'Dark2');

% figure;
% hold on;
% plot(P1, F0, 'o', 'linewidth', 2, 'color', map(3, :));
% xlim([10^-3 10^1]);
% xticks([10^-3 10^-2 10^-1 10^0 10^1]);
% set(gca, 'FontSize', 24, 'XScale', 'log', 'Box', 'on');
% ylabel('Phospholipids syn. flux');
% xlabel('Factor change in rate of ATP use');
% legend({'High LDH, low PDH activity', 'High LDH, high PDH activity', 'Low LDH, high PDH activity'})
% legend boxoff;

hold on;
plot(P1, F1, 'o', 'linewidth', 2, 'color', map(1, :));
xlim([10^-3 10^1]);
xticks([10^-3 10^-2 10^-1 10^0 10^1]);
set(gca, 'FontSize', 16, 'XScale', 'log', 'Box', 'on');
ylabel('ATP use flux');