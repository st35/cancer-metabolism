clc;
clear all;

P1 = logspace(-4, 4, 500);

N_Metabolites = 26;
N_Fluxes = 30;

F0 = zeros(length(P1), 1);
F1 = zeros(length(P1), 1);

Input = zeros(1, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.160 - x0(24); % M_ADP
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
    F1(i) = F(30);
end

F0_b = zeros(length(P1), 1);
F1_b = zeros(length(P1), 1);

Input = zeros(1, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.16*4.0 - x0(24); % M_ADP
x0(26) = 1.0; % C_ATP

for i = length(P1):-1:1
    i
    Input(1) = P1(i);
    [t, x] = ode23tb(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 50000.0], x0);
    x0 = x(end, :)';
    F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

    if sum(isnan(x0)) > 0
        'Dying is easy. Integrating this is hard. - James Wilson.'
    end

    F0_b(i) = F(8);
    F1_b(i) = F(30);
end

map = brewermap(5, 'BrBG');

figure;
hold on;
plot(P1, F0, 'o', 'linewidth', 2, 'color', map(5, :));
plot(P1, F0_b, 'o', 'linewidth', 2, 'color', map(5, :));
xlim([10^-4 10^4]);
xticks([10^-4 10^-2 10^0 10^2 10^4]);
set(gca, 'FontSize', 16, 'XScale', 'log', 'Box', 'on');
ylabel('Phospholipids syn. flux');
xlabel('Factor change in LDH activity');
axes('Position', [.3 .4 .275 .275])
box on;
plot(P1, F1, 'o', 'linewidth', 2, 'color', map(5, :));
hold on;
plot(P1, F1_b, 'o', 'linewidth', 2, 'color', map(5, :));
xlim([10^-4 10^4]);
xticks([10^-4 10^-2 10^0 10^2 10^4])
set(gca, 'FontSize', 16, 'XScale', 'log', 'Box', 'on');
ylabel('ATP:ADP ratio');