
% DDE script - Halyomorpha halys with temperature-dependent delays

% Created by Luca Rossini on 3 November 2023
% Last update 21 July 2025
% e-mail: luca.rossini@ulb.be

% Start to calculate the simulation time

tic

%% Clearing the workspace before the beginning

clear
clc
close all


%% Load the parameters and other inputs

run("Parameters.m")

    % Load the functions from Functions.m

Fun = Functions;

    % Define the error for the model

ModelError = 0.6;

    % Correlation between the stages for the non-diagonal elements of the
    % matrix Q

StageCorrelation = 1;

    % If traps are placed or not and its attraction coefficient

TrapOn = 1;
TrapAttraction = 0.07;

    % Portion of trapped individuals with respect to the total population
    % of the trapped stage within the field.

MortTrap = 0.2;


    %% Map the observations

    % This array is going to indicate which are the states that have
    % been measured. Update H(x) by adding the number of the state that
    % has been observed, x.

    H = [0 0 0 0 0 0 0 0 0 1];

%% Solve the equation - DDE system open loop

    % Define the time span

t_span = load("Parameters.mat", "t_span");
t_span = cell2mat(struct2cell(t_span));


%% Implementation of the Extended Kalman Filter

% Observation error (read from file or vector)
R_errors = ErrAdults;  % Uncertainties from the measurements (real)


%% Time Setup

t0 = t_span(1);
tf = t_span(2);
dt = 1;
tspan = t0:dt:tf;
obs_times = ExpDataDay;  % Takes the observation days from file


%% Initial State (history function)

x0 = InitHist_DDE;


%% Initial EKF values

x_hat = x0; % Initial state
P = eye(10);  % Initial covariance

% Store the solution

X_estimated = zeros(10, length(tspan));
X_estimated(:,1) = x_hat;


%% EKF loop

t_now = 0;

    % Compure the DDE system and store the values to further use into the
    % prediction step below!

solPartial = Fun.dde_solver(t_span, DailyTemp, obs_times, InitHist_DDE, ...
                            SR, FertPar, MortPar_Egg, MortPar_N1, MortPar_N2, ...
                            MortPar_N3, MortPar_N4, MortPar_N5, MortTrap, ...
                            TrapOn, TrapAttraction, DevRate_Egg, DevRate_N1, ...
                            DevRate_N2, DevRate_N3, DevRate_N4, DevRate_N5, ...
                            DevRate_Ad, LagPar_Egg, LagPar_N1, LagPar_N2, ...
                            LagPar_N3, LagPar_N4, LagPar_N5, LagPar_Am, ...
                            LagPar_PreOvi, LagPar_Amf);

% For loop of the EKF

for i = 2:length(tspan)
    t_next = tspan(i);
    dt_i = t_next - t_now;

    % EKF Prediction step - Extraction of the equation value precomputed in
    % solPartial

    x_pred = deval(solPartial, t_next);

    % Injection of the noise on the model, to account the model error and
    % to convert the model in a stochastic model. Q is the process noise
    % matrix.

    Q = Fun.Q_fun(t_now, x_pred, ModelError, StageCorrelation);

    % Compute the Jacobian

    [A, A_lag] = Fun.Jacobian(t_now, x_pred, DailyTemp, obs_times, SR, ...
                                FertPar, MortPar_Egg, MortPar_N1, ...
                                MortPar_N2, MortPar_N3, MortPar_N4, ...
                                MortPar_N5, MortTrap, TrapOn, TrapAttraction,...
                                DevRate_Egg, DevRate_N1, DevRate_N2, ...
                                DevRate_N3, DevRate_N4, DevRate_N5, ...
                                DevRate_Ad);

    % Propagate the covariance

    P = A * P * A' + A_lag * P * A_lag' + Q;

    % Observation update only if in obs_times

    if ismember(t_next, obs_times)
        y = Fun.y_measured(t_next, ExpAdults);
        R = Fun.R_fun(t_next, R_errors(:,1), R_errors(:,2), 10);

        R_scalar = H * R * H';           % 1×1
        S = H * P * H' + R_scalar;       % 1×1
        K = P * H' / S;                  % 10×1

        x_hat = x_pred + K * (y - H * x_pred);
        P = (eye(10) - K * H) * P;
    else
        x_hat = x_pred;
    end

    X_estimated(:,i) = x_hat;
    t_now = t_next;
end


%% Plot the results from 'Plots.m'

run('Plots.m')


%% End the calculation of the simulation time

toc

