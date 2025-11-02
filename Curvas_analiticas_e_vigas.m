%% ========================================================================
%  ANEL PLANO EM ESTADO PLANO DE TENSÕES
%  Comparação: Solução analítica × Teoria simples de viga (energia compl.)
%  Tensões, deformações e deslocamentos
% ========================================================================

clear; clc; close all;

%% ------------------ Dados geométricos e materiais -----------------------

a = 250e-3;            % [m] raio interno
b = 1.2 * a;           % [m] raio externo
E  = 200e9;            % [Pa]
nu = 0.30;             % [-]
t  = 4e-3;             % [m] espessura 
P  = 200;              % [N] força aplicada

Rm = (a + b)/2;        % [m] raio médio
h  = b - a;            % [m] espessura radial
I  = t * h^3 / 12;     % [m^4] momento de inércia (retangular)

% Constantes da função de tensão (item A)
A = 8.2226755804e7;    % [Pa/m³]
B = -4.62525501e5;     % [Pa·m]
D = -2.507916052e7;    % [Pa·m]

% Vetor angular (0 → 2π)
theta = linspace(0, 2*pi, 400);

%% ========================================================================
%  TENSÕES ANALÍTICAS σ_r e σ_θ
% ========================================================================

sigma_r_a = ( 2*A*a - 2*B/a^3 + D/a ) .* cos(theta); % [Pa]
sigma_r_b = ( 2*A*b - 2*B/b^3 + D/b ) .* cos(theta); % [Pa]
sigma_t_a = ( 6*A*a + 2*B/a^3 + D/a ) .* cos(theta); % [Pa]
sigma_t_b = ( 6*A*b + 2*B/b^3 + D/b ) .* cos(theta); % [Pa]

sigma_r_a_MPa = sigma_r_a / 1e6;
sigma_r_b_MPa = sigma_r_b / 1e6;
sigma_t_a_MPa = sigma_t_a / 1e6;
sigma_t_b_MPa = sigma_t_b / 1e6;

%% ------------------ Teoria de viga (corrigida) --------------------------
%  σθ(r=a,θ) = -P*cosθ/(t*h) - 3P(b+a)*cosθ/(t*h^2)
%  σθ(r=b,θ) = -P*cosθ/(t*h) + 3P(b+a)*cosθ/(t*h^2)
sigma_t_beam_a = -P*cos(theta)/(t*h) - 3*P*(b+a)*cos(theta)/(t*h^2);
sigma_t_beam_b = -P*cos(theta)/(t*h) + 3*P*(b+a)*cos(theta)/(t*h^2);
sigma_t_beam_a_MPa = sigma_t_beam_a / 1e6;
sigma_t_beam_b_MPa = sigma_t_beam_b / 1e6;

%% ------------------ Gráfico de tensões --------------------
figure('Color','w','Position',[200 200 850 450]);

% σr (analítica)
subplot(1,2,1)
plot(theta*180/pi, sigma_r_a_MPa, 'b-', 'LineWidth', 1.6); hold on;
plot(theta*180/pi, sigma_r_b_MPa, 'r-', 'LineWidth', 1.6);
xlabel('\theta [graus]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\sigma_r [MPa]', 'FontSize', 12, 'FontWeight', 'bold');
title('Tensão radial \sigma_r (analítica)', 'FontSize', 13);
legend('r=a','r=b','Location','best');
grid on; box on;

% σθ (analítica × teoria de viga)
subplot(1,2,2)
hold on;
% Analítica
plot(theta*180/pi, sigma_t_a_MPa, 'r-', 'LineWidth', 1.6);
plot(theta*180/pi, sigma_t_b_MPa, 'b-', 'LineWidth', 1.6);
% Teoria de viga
plot(theta*180/pi, sigma_t_beam_a_MPa, 'r--', 'LineWidth', 1.4);
plot(theta*180/pi, sigma_t_beam_b_MPa, 'b--', 'LineWidth', 1.4);

xlabel('\theta [graus]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\sigma_\theta [MPa]', 'FontSize', 12, 'FontWeight', 'bold');
title('Tensão circunferencial \sigma_\theta (r=a e r=b)', 'FontSize', 13);
legend({'Analítica r=a','Analítica r=b','Viga r=a','Viga r=b'}, 'Location','best');
grid on; box on; hold off;

sgtitle('Distribuição de Tensões – Analítica × Teoria de Viga (θ)', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
%  DEFORMAÇÕES ANALÍTICAS ε_r e ε_θ
% ========================================================================

eps_r_a = (sigma_r_a - nu .* sigma_t_a) ./ E;
eps_r_b = (sigma_r_b - nu .* sigma_t_b) ./ E;
eps_t_a = (sigma_t_a - nu .* sigma_r_a) ./ E;
eps_t_b = (sigma_t_b - nu .* sigma_r_b) ./ E;

eps_r_a_micro = eps_r_a * 1e6;
eps_r_b_micro = eps_r_b * 1e6;
eps_t_a_micro = eps_t_a * 1e6;
eps_t_b_micro = eps_t_b * 1e6;

% Teoria de viga (εθ = σθ / E)
eps_t_beam_a_micro = (sigma_t_beam_a / E) * 1e6;
eps_t_beam_b_micro = (sigma_t_beam_b / E) * 1e6;

%% ------------------ Gráfico de deformações ----------------
figure('Color','w','Position',[200 200 850 450]);

% εr (analítica)
subplot(1,2,1)
plot(theta*180/pi, eps_r_a_micro, 'b-', 'LineWidth', 1.6); hold on;
plot(theta*180/pi, eps_r_b_micro, 'r-', 'LineWidth', 1.6);
xlabel('\theta [graus]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\epsilon_r [\mue/m]', 'FontSize', 12, 'FontWeight', 'bold');
title('Deformação radial \epsilon_r (analítica)', 'FontSize', 13);
legend('r=a','r=b','Location','best');
grid on; box on;

% εθ (analítica × teoria de viga)
subplot(1,2,2)
hold on;
% Analítica
plot(theta*180/pi, eps_t_a_micro, 'r-', 'LineWidth', 1.6);
plot(theta*180/pi, eps_t_b_micro, 'b-', 'LineWidth', 1.6);
% Teoria de viga
plot(theta*180/pi, eps_t_beam_a_micro, 'r--', 'LineWidth', 1.4);
plot(theta*180/pi, eps_t_beam_b_micro, 'b--', 'LineWidth', 1.4);

xlabel('\theta [graus]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\epsilon_\theta [\mue/m]', 'FontSize', 12, 'FontWeight', 'bold');
title('Deformação circunferencial \epsilon_\theta (r=a e r=b)', 'FontSize', 13);
legend({'Analítica r=a','Analítica r=b','Viga r=a','Viga r=b'}, 'Location','best');
grid on; box on; hold off;

sgtitle('Distribuição de Deformações – Analítica × Teoria de Viga (θ)', ...
    'FontSize', 14, 'FontWeight', 'bold');


%% ========================================================================
%  DESLOCAMENTOS (Energia complementar)
% ========================================================================

% Parâmetro β pequeno 
beta = deg2rad(0.34);

% Fator comum
coef = P*Rm^3*(1 - cos(beta)) / (E*I);

% Deslocamentos pela teoria de viga (energia complementar)
u_theta_viga = coef*((theta - sin(theta)) - (beta - sin(beta)));  % [m]
u_r_viga     = coef*(cos(beta) - cos(theta));                     % [m]
u_theta_viga_mm = u_theta_viga * 1e3;
u_r_viga_mm     = u_r_viga * 1e3;

% Deslocamento analítico (circunferencial, em Rm)
beta_const = -(1/E)*((1 - 3*nu)*A*Rm^2 + (1 + nu)*B/Rm^2 + (1 - nu)*D*log(Rm));
u_theta_analit_mm = (sin(theta)/E) .* ((nu + 5)*A*Rm^2 + (1 + nu)*B/Rm^2 + ...
                      (1 - nu)*(1 - log(Rm))*D) * 1e3 - beta_const*theta*1e3;

%% ------------------ Gráfico de deslocamentos ----------------------------
figure('Color','w','Position',[200 200 950 450]);

subplot(1,2,1)
plot(theta*180/pi, u_r_viga_mm, 'b-', 'LineWidth', 1.6);
xlabel('\theta [graus]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('u_r [mm]', 'FontSize', 12, 'FontWeight', 'bold');
title('Deslocamento radial u_r (teoria de viga)', 'FontSize', 13);
grid on; box on;

subplot(1,2,2)
hold on;
plot(theta*180/pi, u_theta_analit_mm, 'r-', 'LineWidth', 1.6);
plot(theta*180/pi, u_theta_viga_mm, 'k--', 'LineWidth', 1.4);
xlabel('\theta [graus]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('u_\theta [mm]', 'FontSize', 12, 'FontWeight', 'bold');
title('Deslocamento circunferencial u_\theta', 'FontSize', 13);
legend('Analítica','Teoria de viga','Location','best');
grid on; box on; hold off;

sgtitle('Distribuição de Deslocamentos – Analítica × Teoria de Viga (θ)', ...
    'FontSize', 14, 'FontWeight', 'bold');
