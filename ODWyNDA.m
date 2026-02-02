clear all; close all;

% Load data
load NWIData.mat;

% Simulation parameters
K = 3; % Sparsity index (0, 1, 2, or 3)
dt = 0.5; n = 3; 
R = 11; r = R + (K * R);
R1 = 4; R2 = 6; R3 = 1;
r1 = R1 + (K * R1); r2 = R2 + (K * R2); r3 = R3 + (K * R3);

% Model parameters
rho = 1.037; g = 9.81; 
rm = 0.0274; fp = 3; cv = 11.95;
pf1 = 8.45e-2; pf2 = 2.48e-4; pf3 = -7.91e-4; 
pi1 = 0.3075; pi2 = 1.76e-3; pi3 = -4.63e-3;
 
% Variable initialization
xakf = yactArray(:,1); xbar = yactArray(:,1); xhat = yactArray(:,1);
takf = zeros(r,1); tbar = zeros(r,1); that = zeros(r,1);
tb{1} = zeros(r1,1); tb{2} = zeros(r2,1); tb{3} = zeros(r3,1);
theta = [1; rho * g * pf1; rho * g * pf2;  ...
    rho * g * pf3; 1; -rm; -fp; rho * g * pi1; ...
    rho * g * pi2; rho * g * pi3; cv * sqrt(inv(rho))];

% Array initialization
tactArray = [];
xakfArray = []; xbarArray = []; xhatArray = [];
takfArray = []; tbarArray = []; thatArray = [];

% AKF parameters
lambda = 0.999; Q = 0.1*eye(n); R = 0.1*eye(n);
P = 0.01*eye(n); Upsilon = zeros(n,r); S = 1e-8*eye(r); 

% RSR parameters
db = 1e-5; pb = 100;
Vb{1} = pb * eye(r1); Vb{2} = pb * eye(r2); Vb{3} = pb * eye(r3);
Gb{1} = inv(db)*eye(r1); Gb{2} = inv(db)*eye(r2); Gb{3} = inv(db)*eye(r3); 
Ub_prev{1} = inv(Vb{1}); Ub_prev{2} = inv(Vb{2}); Ub_prev{3} = inv(Vb{3});

% RSBL parameters
Sihat = 1e-4 * eye(r); alpha = 0.999;
Ghat_prev = 1e-4 * eye(r); eps = 1e-7;

%% Simulation
for i = 1:length(uactArray)
    % Store variables
    xakfArray = [xakfArray xakf]; takfArray = [takfArray takf];
    xbarArray = [xbarArray xbar]; tbarArray = [tbarArray tbar];
    xhatArray = [xhatArray xhat]; thatArray = [thatArray that];

    % Data indexing
    y = yactArray(:,i); u = uactArray(:,i);

    % Invoke library functions
    [Psi, Phi] = WILibFuncs(K, y, u, pf(i), qf(i));
    tact = [theta(1:R1); zeros(r1 - R1,1); ...
        theta(R1 + 1:R1 + R2); zeros(r2 - R2,1); ...
        theta(R1 + R2 + 1:R1 + R2 + R3); zeros(r3 - R3,1);];
    tactArray = [tactArray tact];
    
    % AKF state covariance update
    P = P + Q;
    Sigma = P + R;
    G = P * inv(Sigma);
    P = (eye(n) - G) * P;

    % AKF parameter covariance update
    Omega = Upsilon + Psi;
    Upsilon = (eye(n) - G) * Upsilon + (eye(n) - G) * Psi;
    Lambda = inv(lambda * Sigma + Omega * S * Omega');
    Gamma = S * Omega' * Lambda;
    S = inv(lambda) * S - inv(lambda) * S * Omega' * Lambda * Omega * S;
    
    % AKF estimation update
    eakf = y - xakf;
    takf = takf + Gamma * eakf;
    xakf = xakf + G * eakf + Upsilon * Gamma * eakf;

    % RSR auxiliary matrix update
    eb = y - xbar;
    for j = 1:n
        for k = 1:length(Phi{j})
            % RSR reweighting penalty
            vb{j}(k) = weib(tb{j}(k));
        end
        Vb{j} = diag(vb{j});
        Ub{j} = inv(pb * Vb{j});
        Gb{j} = Gb{j} - (Gb{j} * Phi{j} * Phi{j}' * Gb{j}) / ...
            (1 + Phi{j}' * Gb{j} * Phi{j});
        Pb{j} = Ub{j} - Ub{j} * inv(Ub{j} + Gb{j}) * Ub{j};
        tb{j} = tb{j} + Pb{j} * Phi{j} * eb(j) - ...
            Pb{j} * (inv(Ub{j}) - inv(Ub_prev{j})) * tb{j};
        Ub_prev{j} = Ub{j};

        % RSR state and parameter update
        xbar(j) = Phi{j}' * tb{j};
    end
    tbar = [tb{1}; tb{2}; tb{3}];

    % RSBL hyperparameter optimization
    ghat = [];
    for j = 1:n
        for k = 1:length(Phi{j})
            gopt = (y(j) * Phi{j}(k) * Phi{j}(k)' * y(j) - ...
                R(j,j) * Phi{j}(k)' * Phi{j}(k)) / ...
                (Phi{j}(k)' * Phi{j}(k) * Phi{j}(k)' * Phi{j}(k) + eps);
            ghat = [ghat; gopt];
        end
    end

    % RSBL parameter update
    Ghat = alpha * Ghat_prev + (1 - alpha) * diag(ghat); 
    DG = Ghat - Ghat_prev; Ghat_prev = Ghat;
    Shat = Psi * Sihat * Psi' + R;
    Khat = Sihat * Psi' * inv(Shat);
    Sihat = Sihat - Khat * Psi * Sihat;

    % RSBL measurement update
    ehat = y - xhat;
    that = that + Khat * ehat - Khat * Psi * Sihat * DG * that;
    xhat = Psi * that;

    % Norm error
    eakf = tact - takf; nakf(i) = norm(eakf);
    ebar = tact - tbar; nbar(i) = norm(ebar); 
    ehat = tact - that; nhat(i) = norm(ehat);
end

%% RMSE calculation
rakf = sqrt(mean(nakf));
rbar = sqrt(mean(nbar));
rhat = sqrt(mean(nhat));

% Print perfromance
fprintf('-------------------------------\n');
fprintf('| Algorithm | Normalized RMSE |\n');
fprintf('-------------------------------\n');
fprintf('| %-9s | %15.7e |\n', 'AKF', rakf);
fprintf('| %-9s | %15.7e |\n', 'RSR', rbar);
fprintf('| %-9s | %15.7e |\n', 'RSBL', rhat);
fprintf('-------------------------------\n');

%% Visualisation
t = [0:0.5:16200] / 60;

% State figure
fh = figure(1);
fh.Position = [150 400 500 400];

subplot(3,1,1);
plot(t, yactArray(1,:), 'k', 'LineWidth', 10); hold on;
plot(t, xakfArray(1,:), '-.', 'LineWidth', 5);
plot(t, xbarArray(1,:), ':', 'LineWidth', 5);
plot(t, xhatArray(1,:), '--', 'LineWidth', 5);
ylabel('$\bar{p_f}\;(\mathrm{bar})$', 'Interpreter','latex');
xlim("tight"); ylim([100 260]); set(gca, 'FontSize', 18); hold off; 

subplot(3,1,2);
plot(t, yactArray(2,:), 'k', 'LineWidth', 10); hold on;
plot(t, xakfArray(2,:), '-.', 'LineWidth', 5);
plot(t, xbarArray(2,:), ':', 'LineWidth', 5);
plot(t, xhatArray(2,:), '--', 'LineWidth', 5);
ylabel('$\bar{p_i}\;(\mathrm{bar})$', 'Interpreter','latex');
xlim("tight"); ylim([200 380]); set(gca, 'FontSize', 18); hold off;  

subplot(3,1,3);
plot(t, yactArray(3,:), 'k', 'LineWidth', 10); hold on;
plot(t, xakfArray(3,:), '-.', 'LineWidth', 5);
plot(t, xbarArray(3,:), ':', 'LineWidth', 5);
plot(t, xhatArray(3,:), '--', 'LineWidth', 5);
legend('Trajectory', 'AKF', 'RSR','RSBL', 'Location', 'best');
xlabel('$t\;(\mathrm{m})$', 'Interpreter','latex');
ylabel('$\bar{q_r}\;(\mathrm{m^3/h})$', 'Interpreter','latex');
xlim("tight"); ylim([-1 460]); set(gca, 'FontSize', 18); hold off; 

%% Function

% Reweighting function
function y = weib(x)
    eb = 0.1; edb = 1e-7;
    y = inv(abs(x) + eb) / sqrt(x^2 + edb);
end

% Sigmoid function
function y = sig(x)
    y = 1 / (1 + exp(-x));
end