clear all; close all;

% RNG choice and load data
rng(100)
load FWIData.mat;

% Simulation parameters
dt = 0.5; 
n = 3; r = 3;

% Model parameters
rho = 1.037; g = 9.81; 
rm = 0.0274; fp = 3; pe = 206.2; cv = 11.95;
hc1 = 8.45e-2; hc2 = 2.48e-4; hc3 = -7.91e-4;
pc1 = -5.95e-3; pc2 = 9.01e-4; pc3 = -1.66e-3;
pi1 = 0.3075; pi2 = 1.76e-3; pi3 = -4.63e-3;

% Variable initialization
theta = rho * g * [hc1; hc2; hc3];
xhat = yactArray(:,1); that = theta;
xhatArray = []; thatArray = [];

% RSBL parameters
Sihat = 1e-4 * eye(r); eta = 0.999; R = 0.1 * eye(n);
Ghat_prev = 1e-4 * eye(r); eps = 1e-7;

%% Simulation
for i = 1:length(uactArray)
    % Store variables
    xhatArray = [xhatArray xhat]; thatArray = [thatArray that];
    
    % Data indexing
    y = yactArray(:,i); u = uactArray(:,i); t = i * dt;

    % Invoke library functions
    Phi = [u(1) * 1e-2; u(1) * qf(i) * 1e-2; qf(i)^2 * 1e-2];
    Psi = [Phi'; zeros(2,r)];

    % RSBL hyperparameter optimization
    ghat = [];
    for k = 1:length(Phi)
        gopt = (y(1) * Phi(k) * Phi(k)' * y(1) - ...
            R(1,1) * Phi(k)' * Phi(k)) / ...
            (Phi(k)' * Phi(k) * Phi(k)' * Phi(k) + eps);
        ghat = [ghat; gopt];
    end

    % RSBL parameter update
    Ghat = eta * Ghat_prev + (1 - eta) * diag(ghat); 
    DG = Ghat - Ghat_prev; Ghat_prev = Ghat;
    Shat = Psi * Sihat * Psi' + R;
    Khat = Sihat * Psi' * inv(Shat);
    Sihat = Sihat - Khat * Psi * Sihat;

    % RSBL excess nonlinear function
    hhat = pump(u(2), qf(i) - xhat(3), pi1, pi2, pi3); 
    fhat = [pf(i); xhat(1) - (rm * qf(i) + fp) + rho * g * hhat * 1e-2; ...
        cv * u(3) * sqrt(max(xhat(1) - pe, 0) * 1e2 / rho)];

    % RSBL measurement update
    ehat = y - xhat;
    that = that + Khat * ehat - Khat * Psi * Sihat * DG * that;
    xhat = fhat + Psi * that;

    % Calculation of degradation index
    [dh(i), dpw(i), deff(i)] = deg_idx(that, uactArray, qf);

    % MLE of RUL model
    if i > 1
        mu(i-1) = mean(diff(deff) / dt);
        si(i-1) = sqrt(mean((diff(deff) - mu(i-1) * dt).^2 / dt));
    end

    % PoF calculation by BNN
    PoF(i) = BBNFuncs(qf(i), u(1), dh(i), dpw(i), deff(i));
end

% RUL prediction
t_obs = [dt:dt:i*dt];
t_sim = [0:dt:80000];
deff_pred = deff;
for k = i:t_sim(end) / dt
    % Wiener process
    deg_pred = mu(end) * k * dt + si(end) * normrnd(0,dt);
    deff_pred = [deff_pred deg_pred];
end

%% Visualisation

% RUL propagation
fh = figure(1);
fh.Position = [500 600 650 300];
plot(t_sim, deff_pred, 'LineWidth', 8); hold on
plot(t_obs, deff, 'k', 'LineWidth', 8);
yline(1, 'r--', {'Failure Threshold'}, 'LabelHorizontalAlignment', ...
    'left', 'LabelVerticalAlignment', 'bottom', 'Linewidth', 5, ...
    'FontSize', 18); xlim("tight"), ylim("tight");
xline(t_sim(min(find(deff_pred > 1))), 'k-.', ...
    {'Predicted Failure Time'}, 'LabelVerticalAlignment', 'bottom', ...
    'LineWidth', 5, 'FontSize', 18); grid("on");
legend('Prediction', 'Observation', 'Location', 'best')
xlabel('$t\;(\mathrm{h})$', 'Interpreter','latex');
ylabel('$d$', 'Interpreter','latex');  
set(gca, 'FontSize', 18); hold off

% PDF evolution
fh = figure(2);
fh.Position =[600 500 550 350];
t = [50000:1000:550000]; 

% PDF calculation
obsp = 5000; obsp0 = obsp; obsd = 20;
for k = 1:obsd
    % Wiener PDF
    pdf(k,:) = (sqrt(2 * pi * si(obsp).^2 .* t.^3)).^(-1) .* ...
         exp(-(1 - mu(obsp) .* t).^2 ./ (2 * si(obsp).^2 .* t));
    obs_pdf = obsp * ones(size(t)); 

    % Plot PDF
    h1 = plot3(t, obs_pdf , pdf(k,:), 'k', 'LineWidth', 3); 
    hold on;
    
    % MTTF calculation
    [max_pdf(k), idx(k)] = max(pdf(k,:));
    max_pred(k) = t(idx(k)); max_obs(k) = obsp;

    % Update observation period
    obsp = round(obsp + (i - obsp0 * dt) / obsd);
end

% Plot MTTF
h2 = plot3(max_pred, max_obs, max_pdf, 'o-', ...
    'LineWidth', 5, 'MarkerSize', 10); grid("on");
legend([h2], {'MTTF'}, 'Location', 'best');
xlim("tight"), ylim("tight"); zlim("tight");
xlabel('$\mathrm{Simulation\;Time\;(h)}$', 'Interpreter','latex');
ylabel('$\mathrm{Observation\;Period\;(h)}$', 'Interpreter','latex'); 
zlabel('$\mathrm{PDF}$', 'Interpreter','latex');
view([45, 30]); set(gca, 'FontSize', 18); hold off;

% PoF plot by BBN
fh = figure(3);
fh.Position =[600 500 650 300];
plot(t_obs, PoF, 'LineWidth', 10);
xlabel('$t\;(\mathrm{h})$', 'Interpreter','latex');
ylabel('$P(F=1)$', 'Interpreter','latex');
grid on; xlim tight, ylim tight;
set(gca, 'FontSize', 18); hold off;

%% Function

% Degradation evidence generation function
function [dh, dpw, deff] = deg_idx(theta, u, q)
    rho = 1.037; g = 9.81; 
    dh_max = 0.15; dpw_max = 0.15; deff_max = 0.2;
    hc1 = 8.45e-2; hc2 = 2.48e-4; hc3 = -7.91e-4;
    pc1 = -5.95e-3; pc2 = 9.01e-4; pc3 = -1.66e-3;

    u0 = u(1,1); q0 = q(1);
    h0 = pump(u0, q0, hc1, hc2, hc3); 
    pw0 = rho * g * q0 * h0 / 3600;
    ps = pump(u0, q0, pc1, pc2, pc3);
    e0 = pw0 / ps; 

    h = pump(u0, q0, theta(1), theta(2), theta(3));
    pw = q0 * h / 3600;
    eff = pw / ps;

    dh = (h0 - h / (rho * g)) / (h0 * dh_max);
    dpw = (pw0 - pw) / (pw0 * dpw_max);
    deff = (e0 - eff) / deff_max;
end

% Pump curve function
function hp = pump(u, q, p1, p2, p3)
    hp = p1 * u + p2 * u * q + p3 * q^2;
end