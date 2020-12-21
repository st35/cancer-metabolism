clc;
clear all;

P1 = logspace(-3, 0, 500);

N_Metabolites = 31;
N_Fluxes = 38;

F0 = zeros(length(P1), 1);
S0 = zeros(length(P1), 1);
S1 = zeros(length(P1), 1);
S2 = zeros(length(P1), 1);
S3 = zeros(length(P1), 1);

Input = zeros(1, 2);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(22) = 1.0e-2; % M_NADH
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.16*4.0 - x0(24); % M_ADP
x0(26) = 1.0; % C_ATP

Input(2) = 1.0e-4;

for i = length(P1):-1:1
    i
    Input(1) = P1(i);
    [t, x] = ode23tb(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 50000.0], x0);
    x0 = x(end, :)';
    F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

    if sum(isnan(x0)) > 0
        'Dying is easy. Integrating this is hard. - James Wilson.'
    end

    F0(i) = F(38);

    S0(i) = x0(16);
    S1(i) = x0(18);
    S2(i) = x0(21);
    S3(i) = (0.3 - x0(22)) / x0(22);
end

nflux = F0;
nconc = [S0(500) S1(500) S2(500) S3(500)]; % P1(167) = 0.01

F0 = zeros(length(P1), 1);
S0 = zeros(length(P1), 1);
S1 = zeros(length(P1), 1);
S2 = zeros(length(P1), 1);
S3 = zeros(length(P1), 1);

Input = zeros(1, 2);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(22) = 1.0e-2; % M_NADH
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.16*4.0 - x0(24); % M_ADP
x0(26) = 1.0; % C_ATP

Input(2) = 1.0e4;

for i = length(P1):-1:1
    i
    Input(1) = P1(i);
    [t, x] = ode23tb(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 50000.0], x0);
    x0 = x(end, :)';
    F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

    if sum(isnan(x0)) > 0
        'Dying is easy. Integrating this is hard. - James Wilson.'
    end

    F0(i) = F(38);

    S0(i) = x0(16);
    S1(i) = x0(18);
    S2(i) = x0(21);
    S3(i) = (0.3 - x0(22)) / x0(22);
end

rflux = F0;
rconc = [S0(167) S1(167) S2(167) S3(167)]; % P1(167) = 0.01

map = brewermap(3, 'Dark2');
figure;
hold on;
plot(P1, nflux, 'o', 'linewidth', 2, 'color', map(1, :));
plot(P1, rflux, 'o', 'linewidth', 2, 'color', map(2, :));
set(gca, 'FontSize', 16, 'Box', 'on');
set(gca, 'XScale', 'log');
legend({'Low NADPH-dependent IDH activity', 'High NADPH-dependent IDH activity'});
legend boxoff;
legend('location', 'northwest');
ylabel('Fatty acids syn. flux');
xlabel('O_{2} conc.');

map = brewermap(4, 'Paired');
figure;
b = bar(1:4, log2(nconc ./ rconc), 'EdgeAlpha', 0, 'FaceColor', 'flat');
for i = 1:4
    b.CData(i, :) = map(i, :);
end
set(gca, 'FontSize', 16, 'Box', 'on');
xticks([]);
ylabel(['Fold change (log_{2}) in the' 10 'mitochondrial conc.']);