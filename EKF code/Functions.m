
% Definition of the functions contained in the RunmeDDE_HalyomorphaVD.m 
% file

% Created by Luca Rossini on 3 November 2023
% Last update 21 July 2025
% E-mail: luca.rossini@ulb.be


classdef Functions

    methods (Static) % Insert any function in between "methods" and "end"


    %% Temperature reader

        function Temp = TempFunction(time, TempArray)
            
            % This if/else statement is needed because there is a problem
            % when time is 0, because it goes to read the array in the 0
            % position while Matlab starts to count from 1

            if time <= 0
                Temp = TempArray(1);
            else
                Temp = TempArray(time);
            end 
        end


    %% Delay function 1 - Rude interpolation of the minimum delays case of
    % the Egg and N1

    function Delay = DelayFun_EggN1(T, Parameters)

        Delay = (Parameters(1) * T^2 + Parameters(2) * T + ...
                   Parameters(3)) / (T + Parameters(4));

        if Delay <= 0
            Delay = 1;
            
        elseif T <= 15 % If out of the lower threshold

            T = 15;
            Delay = (Parameters(1) * T^2 + Parameters(2) * T + ...
                        Parameters(3)) / (T + Parameters(4));

        else
            Delay = (Parameters(1) * T^2 + Parameters(2) * T + ...
                   Parameters(3)) / (T + Parameters(4));
        end
        
    end


    %% Delay function 2 - Rude interpolation of the minimum delays case of
    % the N2

    function Delay = DelayFun_N2(T, Parameters)

        Delay = Parameters(1) * exp(Parameters(2) * T);

        if Delay <= 0
            Delay = 1;

        elseif T <= 15 % If out of the lower threshold

            T = 15;
            Delay = Parameters(1) * exp(Parameters(2) * T);

        else
            Delay = Parameters(1) * exp(Parameters(2) * T);
        end
    end


    %% Delay function 3 - Rude interpolation of the minimum delays case of
    %  the N3

    function Delay = DelayFun_N3(T, Parameters)

        Delay = Parameters(1) * T^2 + Parameters(2) * T + Parameters(3);

        if Delay <= 0
            Delay = 1;
        else
            Delay = Parameters(1) * T^2 + Parameters(2) * T ...
                    + Parameters(3);
        end
    end


    %% Delay function 4 - Rude interpolation of the minimum delays case of
    %  the N4

    function Delay = DelayFun_N4(T, Parameters)

        Delay = Parameters(1) * T^3 + Parameters(2) * T^2 + ...
                 Parameters(3) * T + Parameters(4);

        if Delay <= 0
            Delay = 1;
        else
            Delay = Parameters(1) * T^3 + Parameters(2) * T^2 + ...
                     Parameters(3) * T + Parameters(4);
        end
    end


    %% Delay function 5 - Rude interpolation of the minimum delays case of
    %  the N5

    function Delay = DelayFun_N5(T, Parameters)

        Delay = (Parameters(1) * T + Parameters(2)) / (T + Parameters(3));

        if Delay <= 0
            Delay = 1;
        else
            Delay = (Parameters(1) * T + Parameters(2)) / ...
                     (T + Parameters(3));
        end
    end


    %% Delay function 6 - Rude interpolation of the preoviposition period

    function Delay = DelayFun_PreOvi(T, Parameters)

        Delay = Parameters(1) * T^Parameters(2) + Parameters(3);

        if Delay <= 0
            Delay = 1;
        else
            Delay = Parameters(1) * T^Parameters(2) + Parameters(3);
        end
    end


    %% Temperature-dependent delays - To feed the DDE

    function Del = Delays(t, ~, TempArray, LagPar_Egg, LagPar_N1, ...
                          LagPar_N2, LagPar_N3, LagPar_N4, LagPar_N5, ...
                          LagPar_Am, LagPar_PreOvi, LagPar_Amf)

            % Import the parameters and functions

            F = Functions;

            % Manage time and temperatures

            time = int32(t);
            DayTemp = F.TempFunction(time, TempArray);

            tau_e = t - F.DelayFun_EggN1(DayTemp, LagPar_Egg);  % Eggs lag
            tau_N1 = t - F.DelayFun_EggN1(DayTemp, LagPar_N1); % N1 lag
            tau_N2 = t - F.DelayFun_N2(DayTemp, LagPar_N2); % N2 lag
            tau_N3 = t - F.DelayFun_N3(DayTemp, LagPar_N3); % N3 lag
            tau_N4 = t - F.DelayFun_N4(DayTemp, LagPar_N4); % N4 lag
            tau_N5 = t - F.DelayFun_N5(DayTemp, LagPar_N5); % N5 lag
            tau_Am = t - LagPar_Am(1);       % Males lag - There is no lag
            tau_Anmf = t - F.DelayFun_PreOvi(DayTemp, LagPar_PreOvi);   
                                            % Female mating lag
            tau_Amf = t - LagPar_Amf(1);     % Female lag - There is no lag

            tau_Trap = 1;     % Trapped individuals

            Del = [tau_e, tau_N1, tau_N2, tau_N3, tau_N4, tau_N5, ...
                   tau_Am, tau_Anmf, tau_Amf, tau_Trap];
             
        end


    %% Fertility rate function

        function B = FertRate(T, Parameters)

            alpha = Parameters(1); 
            gamma = Parameters(2); 
            lambda = Parameters(3); 
            delta = Parameters(4); 
            tau = Parameters(5);
            
            B = alpha * ((gamma + 1)/(pi * lambda^(2 * gamma + 2)) * ...
                   (lambda^2 - ((T - tau)^2 + delta^2))^gamma);

        end


    %% Development rate function

        function G = DevRate(T, Parameters)

            a = Parameters(1);
            T_L = Parameters(2);
            T_M = Parameters(3);
            m = Parameters(4);

            % This if/else statement is needed to ensure that the solution
            % is biologically reasonable
            
            if T < T_L || T > T_M
                G = 0;
            else
                G = a * T * (T - T_L) * (T_M - T)^(1/m);
            end

        end


    %% Mortality rate function

        function M = MortRate(T, Parameters)

            a = Parameters(1); 
            b = Parameters(2);
            c = Parameters(3);
            d = Parameters(4);
            e = Parameters(5);

            M = a * T^4 + b * T^3 + c * T^2 + d * T + e;

            % This if/else statement is needed to ensure that the solution
            % is biologically reasonable
            
            if M < 0
                M = 0;
            else
                M = (a * T^4 + b * T^3 + c * T^2 + d * T + e);
            end
               
        end
       

    %% DDE solver wrapper

    % This function takes the parameters on the DDE as an input and
    % provides the solution of the system as an output. This is helpful to
    % integrate the DDE over custom ranges, avoiding the implementation by
    % hand

        function sol = dde_solver(t_span, DailyTemp, ObsTimes, InitHist_DDE, ...
                                  SR, FertPar, MortPar_Egg, MortPar_N1, ...
                                  MortPar_N2, MortPar_N3, MortPar_N4, ...
                                  MortPar_N5, MortTrap, TrapOn, TrapAttraction, ...
                                  DevRate_Egg, DevRate_N1, ...
                                  DevRate_N2, DevRate_N3, DevRate_N4, ...
                                  DevRate_N5, DevRate_Ad, LagPar_Egg, ...
                                  LagPar_N1, LagPar_N2, LagPar_N3, ...
                                  LagPar_N4, LagPar_N5, LagPar_Am, ...
                                  LagPar_PreOvi, LagPar_Amf)
            
            % Calling the functions contained into the static class
            % Functions

            F = Functions;
            
            % Solving the DDE system

            sol = ddesd(@(t,y,Z) F.ddefun_partial(t,y,Z,SR,FertPar, ...
                                MortPar_Egg, MortPar_N1, MortPar_N2, ...
                                MortPar_N3, MortPar_N4, MortPar_N5, MortTrap, ...
                                TrapOn, TrapAttraction,...
                                DevRate_Egg, DevRate_N1, DevRate_N2, ...
                                DevRate_N3, DevRate_N4, DevRate_N5, ...
                                DevRate_Ad, DailyTemp, ObsTimes), ...
                        @(t,y) F.Delays(t,y,DailyTemp, LagPar_Egg, ...
                                LagPar_N1, LagPar_N2,LagPar_N3, LagPar_N4, ...
                                LagPar_N5, LagPar_Am, LagPar_PreOvi, ...
                                LagPar_Amf), ...
                        @(t) F.history(t, InitHist_DDE), t_span);

        end


    %% DDE system 

    % This is the declaration of the DDE system of equations

        function dydt_partial = ddefun_partial(t, y, Z, SR, FertPar, ...
                                MortPar_Egg, MortPar_N1, MortPar_N2, ...
                                MortPar_N3, MortPar_N4, MortPar_N5, ...
                                MortTrap, TrapOn, TrapAttraction, DevRate_Egg, DevRate_N1, ...
                                DevRate_N2, DevRate_N3, DevRate_N4, ...
                                DevRate_N5, DevRate_Ad, TempArray, ObsTimes)

            % Approximation of the lags - Needed by ddesd even for the
            % states were the delay is not foreseen!

            ylag1 = Z(1);    % Eggs lag
            ylag2 = Z(2);    % N1 lag
            ylag3 = Z(3);    % N2 lag
            ylag4 = Z(4);    % N3 lag
            ylag5 = Z(5);    % N4 lag
            ylag6 = Z(6);    % N5 lag
            ylag7 = Z(7);    % Adult males lag
            ylag8 = Z(8);    % Female mating lag - Preoviposition period!
            ylag9 = Z(9);    % Female lag
            ylag10 = Z(10);  % Trapped adults

            % It needs to call the class containing the functions needed,
            % also if the present function is contained in the same class!

            F = Functions;

            % Calculate the daily temperature from the
            % 'TemperatureInput.xlsx' file in 'Parameters.m'

            time = int32(t);
            temp = F.TempFunction(time, TempArray);

            % Check the observation times: empties the trap

            if any(ObsTimes == time)
                R_T = 1;%TrapAttraction;
            else
                R_T = 0;
            end

            % DDE system
        
                % y(1) = eggs without lag
                % y(2) = N1 without lag
                % y(3) = N2 without lag
                % y(4) = N3 without lag
                % y(5) = N4 without lag
                % y(6) = N5 without lag                
                % y(7) = males without lag
                % y(8) = non mated females without lag
                % y(9) = mated females without lag
                % y(10) = Trapped adults without lag
            
            dydt_partial = [ F.DevRate(temp, DevRate_Ad) * ...             % Eggs
                               F.FertRate(temp, FertPar) * y(9) ...
                             - F.DevRate(temp, DevRate_Egg) * y(1) ...
                             - F.MortRate(temp, MortPar_Egg) * y(1);

                             F.DevRate(temp, DevRate_Egg) * ylag1 ...      % N1
                             - F.DevRate(temp, DevRate_N1) * y(2) ...
                             - F.MortRate(temp, MortPar_N1) * y(2);

                             F.DevRate(temp, DevRate_N1) * ylag2 ...       % N2
                             - F.DevRate(temp, DevRate_N2) * y(3) ...
                             - F.MortRate(temp, MortPar_N2) * y(3);

                             F.DevRate(temp, DevRate_N2) * ylag3 ...       % N3
                             - F.DevRate(temp, DevRate_N3) * y(4) ...
                             - F.MortRate(temp, MortPar_N3) * y(4);

                             F.DevRate(temp, DevRate_N3) * ylag4 ...       % N4
                             - F.DevRate(temp, DevRate_N4) * y(5) ...
                             - F.MortRate(temp, MortPar_N4) * y(5);

                             F.DevRate(temp, DevRate_N4) * ylag5 ...       % N5
                             - F.DevRate(temp, DevRate_N5) * y(6) ...
                             - F.MortRate(temp, MortPar_N5) * y(6);

                             (1 - SR) * F.DevRate(temp, DevRate_N5) * ...  % Am 
                             ylag6 - F.DevRate(temp, DevRate_Ad) * y(7)...
                             - TrapOn * TrapAttraction * y(7);

                             SR * F.DevRate(temp, DevRate_N5) * ylag6 ...  % Anmf
                             - y(8) - TrapOn * TrapAttraction * y(8);

                             y(8) - F.DevRate(temp, DevRate_Ad) * ...     % Amf
                             y(8) - F.DevRate(temp, DevRate_Ad) * y(9) ...
                             - TrapOn * TrapAttraction * y(9);
                             
                             TrapOn * y(7) * MortTrap + ...               % Trapped individuals
                             TrapOn * y(8) * MortTrap + ...
                             TrapOn * y(9) * MortTrap ...
                             - R_T * y(10)];

        end


    %% Define the initial history for the DDE

    % This is the equivalent of the initial conditions for an ODE system

        function s = history(t, s)

            s = 2 * heaviside(t) * s;

        end


    %% Compute the Jacobian matrix

        function [J_now, J_lags] = Jacobian(t, y, DailyTemp, ObsTimes, SR, ...
                                    FertPar, MortPar_Egg, MortPar_N1, ...
                                    MortPar_N2, MortPar_N3, MortPar_N4, ...
                                    MortPar_N5, MortTrap, TrapOn, TrapAttraction,...
                                    DevRate_Egg, DevRate_N1, DevRate_N2, ...
                                    DevRate_N3, DevRate_N4, DevRate_N5, ...
                                    DevRate_Ad)

            % Define the functions

            F = Functions;
            temp = F.TempFunction(int32(t), DailyTemp);

            % Account the trap inspections

            if any(ObsTimes == t)
                R_T = 1;%TrapAttraction;
            else
                R_T = 0;
            end

            % Define the matrices that store the values of the Jacobian at
            % time t and at time t-tau

            J_now = zeros(10, 10); % Current time Jacobian
            J_lags = zeros(10, 10); % Delayed Jacobian

            % Jacobian with respect to y(t)
                
                % Egg and Fertility - Equation 1
            J_now(1,1) = -F.DevRate(temp, DevRate_Egg) ...
                         - F.MortRate(temp, MortPar_Egg);
            J_now(1,9) = F.FertRate(temp, FertPar) * F.DevRate(temp, DevRate_Ad);

                % N1 - Equation 2
            J_now(2,2) = -F.DevRate(temp, DevRate_N1) ...
                         - F.MortRate(temp, MortPar_N1);

                % N2 - Equation 3
            J_now(3,3) = -F.DevRate(temp, DevRate_N2) ...
                         - F.MortRate(temp, MortPar_N2);

                % N3 - Equation 4
            J_now(4,4) = -F.DevRate(temp, DevRate_N3) ...
                         - F.MortRate(temp, MortPar_N3);

                % N4 - Equation 5
            J_now(5,5) = -F.DevRate(temp, DevRate_N4) ...
                         - F.MortRate(temp, MortPar_N4);

                % N5 - Equation 6
            J_now(6,6) = -F.DevRate(temp, DevRate_N5) ...
                         - F.MortRate(temp, MortPar_N5);

                % Am - Equation 7
            J_now(7,7) = -F.DevRate(temp, DevRate_Ad) - TrapOn * TrapAttraction;          
            
                % Anmf - Equation 8
            J_now(8,8) = -1 - TrapOn * TrapAttraction;                  

                % Amf - Equation 9
            J_now(9,8) = 1 - F.DevRate(temp, DevRate_Ad);

            J_now(9,9) = -F.DevRate(temp, DevRate_Ad) - TrapOn * TrapAttraction;

                % Trapped males
            J_now(10,7) = TrapOn * MortTrap;

                % Trapped non mated females
            J_now(10,8) = TrapOn * MortTrap;

                % Trapped mated females
            J_now(10,9) = TrapOn * MortTrap;

                % Trap emptied
            J_now(10,10) = -R_T;

            % Jacobian with respect to y(t - τi)

            % Rows correspond to the equation, columns correspond to the 
            % delayed variable

            J_lags(2, 1) = F.DevRate(temp, DevRate_Egg);  % ylag1 = Egg
            J_lags(3, 2) = F.DevRate(temp, DevRate_N1);   % ylag2 = N1
            J_lags(4, 3) = F.DevRate(temp, DevRate_N2);   % ylag3 = N2
            J_lags(5, 4) = F.DevRate(temp, DevRate_N3);   % ylag4 = N3
            J_lags(6, 5) = F.DevRate(temp, DevRate_N4);   % ylag5 = N4

            J_lags(7, 6) = (1 - SR) * F.DevRate(temp, DevRate_N5); % ylag6 = Am
            J_lags(8, 6) = SR * F.DevRate(temp, DevRate_N5); % ylag7 = Anmf

        end 


    %% Process noise matrix Q(x)

    % Uncertainty associated with the model - RANDOM NOISE OF THE MODEL!

    function Q = Q_fun(t, x, ErrPerc, CorrelationPercentage)
    
            % Set epsilon to avoid zeroes

            epsilon = 1e-6;
    
            % Variance proportional to the population into the compartment
            
            q_diag = (ErrPerc * abs(x) + epsilon).^2;
    
            Q = diag(q_diag);

                % Insert the correlation between the trapped stage (the one
                % observed) and the other stages. If the matrix Q is
                % diagonal (the last line of code above) we are assuming no
                % correlation between the stages, and the correction on one
                % state is not always correcting the others properly.

            sigmaT   = sqrt(Q(10,10));

            for i = 1:9
                sigma_i    = sqrt(Q(i,i));
                Q(i,10)    = CorrelationPercentage * sigma_i * sigmaT;
                Q(10,i)    = Q(i,10);
            end

        end


    %% Measurement noise matrix R(t)

    % Uncertainty associated with the measurements

    function R = R_fun(t, obs_times, obs_error_var, n_state)
        
        % R_fun: provides the covariance matrix of the observations at time
        % t
        %
        % INPUTS:
        %   t               - current time
        %   obs_times       - array of the observations
        %   obs_error_var   - array of the variance/uncertainties
        %                     associated with the observations
        %   n_state         - number of state variables (dimension of R)
        %
        % OUTPUT:
        %   R               - diagonal matrix n_state x n_state, with  
        %                     error on the states observed

        % Interpolation of the variance over time

        var_t = interp1(obs_times, obs_error_var, t, 'linear', 'extrap');

        % Building the matrix R: Diagonal with zeroes, except on the states
        % where measurements are available

        R = zeros(n_state);
        R(10,10) = var_t;
    end


    %% EKF Intermittent state update

    % This function is going to update the estimation when a new
    % measurement is available from the field

    function [x_new, P_new] = EKF_update_dde(x_pred, P_pred, y_meas, t_k, ...
                                             R_fun, H_fun, t_meas)
    
        % Controlla se il tempo corrente t_k è tra i tempi di misura
        % Check if the current time t_k is available between the
        % measurement times

        tol = 1e-6;  % tolerance threshold for time comparison
    
        if any(abs(t_meas - t_k) < tol)

            % If there is a new measurement, update the estimation

            Hk = H_fun(x_pred, t_k);
            Rk = R_fun(t_k);
            y_hat = Hk * x_pred;

            % Residual

            innovation = y_meas - y_hat;

            % Kalman gain

            S = Hk * P_pred * Hk' + Rk;
            K = P_pred * Hk' / S;

            % Updated state

            x_new = x_pred + K * innovation;

            % Updated covariance

            P_new = (eye(size(P_pred)) - K * Hk) * P_pred;
        else
            % If there are no measurements, skip the update
            x_new = x_pred;
            P_new = P_pred;
        end

    end


     %% Provides the field measurement at time t

     function y = y_measured(t, Measurements)

        % Finds the value in ExpAdults corresponding to t

        idx = find(Measurements(:,1)==t,1);
        if isempty(idx)
            y = NaN;
        else
            y = Measurements(idx,2);
        end
     end

     
    end

end
