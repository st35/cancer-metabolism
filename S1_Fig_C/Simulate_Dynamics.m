clc;
clear all;

P1 = logspace(-2, 2, 500);

N_Metabolites = 28;
N_Fluxes = 32;

F0 = zeros(length(P1), 1);

Fall = zeros(2, length(P1));

options = odeset('AbsTol', 1.0e-9);

Input = zeros(1, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.160 - x0(24); % M_ADP
x0(26) = 1.0; % C_ATP
x0(28) = 0.001; % C_NADH

for i = 1:length(P1)
    i
    Input(1) = P1(i); % ATP use
    Input(2) = 5.0e2; % LDH
    Input(3) = 1.0; % PDH

    [t, x] = ode23tb(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 500000.0], x0, options);
    x0 = x(end, :)';
    F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

    if sum(isnan(x0)) > 0
        'Dying is easy. Integrating this is hard. - James Wilson.'
    end

    F0(i) = F(30); % NAD use
    Fall(1, i) = F0(i);
end

Input = zeros(1, 1);
x0 = zeros(N_Metabolites, 1);

x0(21) = 4.0e-2*1e2; % M_OAA
x0(24) = 900.0e-4; % M_ATP
x0(25) = 4.160 - x0(24); % M_ADP
x0(26) = 1.0; % C_ATP
x0(28) = 0.001; % C_NADH

for i = 1:length(P1)
    i
    Input(1) = P1(i); % ATP use
    Input(2) = 1.0; % LDH
    Input(3) = 5.0e2; % PDH

    [t, x] = ode23tb(@(t, x) Metabolic_System(0, N_Metabolites, N_Fluxes, Input, x), [0 500000.0], x0, options);
    x0 = x(end, :)';
    F = Metabolic_System(1, N_Metabolites, N_Fluxes, Input, x0)';

    if sum(isnan(x0)) > 0
        'Dying is easy. Integrating this is hard. - James Wilson.'
    end

    F0(i) = F(30); % NAD use
    Fall(2, i) = F0(i);
end

figure;
hold on;
plot(P1, Fall(1, :), 'bo', 'linewidth', 2);
plot(P1, Fall(2, :), 'ro', 'linewidth', 2);
set(gca, 'FontSize', 16, 'XScale', 'log', 'Box', 'on');
xlabel('Factor change in rate of ATP use');
ylabel('NAD use flux');
xticks(logspace(-2, 2, 5));
legend({['High LDH +' 10 'low PDH'], ['Low LDH +' 10 'high PDH']}, 'location', 'northwest');