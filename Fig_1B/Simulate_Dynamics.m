clc;
clear all;

P1 = logspace(-1, 3, 500);

N_Metabolites = 21;
N_Fluxes = 24;

F0 = zeros(length(P1), 1);
F1 = zeros(length(P1), 1);

Input = zeros(1, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA

for i = length(P1):-1:1
    i
    Input(1) = P1(i);
    [t, x] = ode23s(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 500000.0], x0);
    x0 = x(end, :)';
    F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

    if sum(isnan(x0)) > 0
        'Dying is easy. Integrating this is hard. - James Wilson.'
    end

    F0(i) = F(8);
    F1(i) = F(4);
end

map = brewermap(4, 'Paired');

fig = figure;
set(fig, 'defaultAxesColorOrder', [[0 0 0]; [0 0 0]]);
hold on;
yyaxis left;
h1 = plot(P1, F0, 'o', 'linewidth', 2, 'color', map(2, :));
set(gca, 'FontSize', 16, 'XScale', 'log', 'Box', 'on');
xlabel('Factor change in rate of glucose import');
ylabel('Phospholipids syn. flux');
yyaxis right;
h2 = plot(P1, F1, 'o', 'linewidth', 2, 'color', map(4, :));
set(gca, 'FontSize', 16, 'XScale', 'log', 'Box', 'on');
ylabel('Ribose syn. flux');
legend([h1 h2], {'Phospholipids syn. flux', 'Ribose syn. flux'});
legend('location', 'northwest');
legend boxoff;
xticks(logspace(-1, 3, 5));