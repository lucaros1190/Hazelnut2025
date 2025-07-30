
% File that plots the results calculated in 'RunmeDDE_HalyomorphaVD.m'

% Created by Luca Rossini on 27 November 2023
% Last update 5 March 2024
% e-mail: luca.rossini@unitus.it


%% Compute the total adults - Open loop

TotAdultsDDE_OpenLoop = solPartial.y(7, :) + solPartial.y(8, :) ...
                        + solPartial.y(9, :);

TotNymphsDDE_OpenLoop = solPartial.y(6, :);


%% Compute the total adults - EKF

TotAdultsDDE_EKF = X_estimated(7, :) + X_estimated(8, :) ...
                        + X_estimated(9, :);

TotNymphsDDE_EKF = X_estimated(6, :);


%% Plot the solutions of adults and nymphs vs experimental data

figure

subplot(2, 1, 1)

hold on

plot(solPartial.x, TotNymphsDDE_OpenLoop, '--', 'LineWidth', 1.5, 'Color', ...
     [0.4940 0.1840 0.5560])
plot(tspan, TotNymphsDDE_EKF, 'LineWidth', 1.5, 'Color', ...
     [0.4660 0.6740 0.1880])
scatter(ExpDataDay, ExpDataNymphs, 'Marker', '*', 'MarkerEdgeColor', ...
        'black', 'MarkerFaceColor', 'black')
title('DDE model - Nymphs')
xlabel('Time (days)')
ylabel('Population density')
legend('Est nymphs - Open Loop', 'Est nymphs - EKF', 'Exp nymphs')

hold off


subplot(2, 1, 2)

hold on

plot(solPartial.x, TotAdultsDDE_OpenLoop, '--', 'LineWidth', 1.5, 'Color', ...
     [0.4940 0.1840 0.5560])
plot(tspan, TotAdultsDDE_EKF, 'LineWidth', 1.5, 'Color', ...
     [0.4660 0.6740 0.1880])
scatter(ExpDataDay, ExpDataAdults, 'Marker', '*', 'MarkerEdgeColor', ...
        'black', 'MarkerFaceColor', 'black')
title('DDE model - Adults')
xlabel('Time (days)')
ylabel('Population density')
legend('Est adults - Open loop', 'Est adults - EKF', 'Exp adults')

hold off


%% Plot the trapped individuals

figure

hold on

plot(solPartial.x, solPartial.y(10,:), '--', 'LineWidth', 1.5, 'Color', ...
     [0.4940 0.1840 0.5560])
plot(tspan, X_estimated(10,:),  'LineWidth', 1.5, 'Color', ...
     [0.4660 0.6740 0.1880])
scatter(ExpDataDay, ExpDataAdults, 'Marker', '*', 'MarkerFaceColor', ...
        'black')
title('Trapped individuals')
xlabel('Time (days)')
ylabel('Adult trapped')
legend('Estimated- Open loop', 'Estimated - EKF', 'Trapped')