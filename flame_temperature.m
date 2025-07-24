function T_p = flame_temperature(m, n, hf_o, T_air, exc_air, RH, humidity, Q_in, W_in, Q_out, W_out)
% flame_temperature  Computes the final temperature of combustion products (T_p)
%
% INPUTS (in order):
%   m        - Number of carbon atoms (C) in the fuel
%   n        - Number of hydrogen atoms (H) in the fuel
%   hf_o     - Enthalpy of formation of the fuel (kJ/mol)
%   T_air    - Inlet air temperature (°C)
%   exc_air  - Air excess ratio
%   RH       - Relative humidity (between 0 and 1)
%   humidity - Boolean: true if humidity is considered, false otherwise
%   Q_in     - Heat input (kJ/kmol)
%   W_in     - Work input (kJ/kmol)
%   Q_out    - Heat output (kJ/kmol)
%   W_out    - Work output (kJ/kmol)
%
% OUTPUT:
%   T_p  - Solution for the final temperature of combustion products
%          or adiabatic flame temperature (K)
%   Mário Loureiro - 2025/ a2021112398@isec.pt

    syms T_p

    T_airK = T_air + 273.15;
    p_atm = 101.325e3; % Pa

    a = (m + n/4) * exc_air;
    b = (m + n/4) * 3.76 * exc_air;

    if humidity
        p_v = RH * XSteam('psat_T', T_air) * 1e5;
        c = ((p_v / p_atm) * ((m + n/4) * 4.76 * exc_air)) / (1 - (p_v / p_atm));
    else
        c = 0;
    end

    d = m;
    e = n/2 + c;
    f = b;
    g = a - (m + n/4);

    % Reactants enthalpies at T_airK
    h_o2  = -9e-17*T_airK^6 + 9e-13*T_airK^5 - 3e-9*T_airK^4 + 4e-6*T_airK^3 + 0.0019*T_airK^2 + 27.751*T_airK + 126.42;
    h_n2  = 3e-12*T_airK^4 - 5e-7*T_airK^3 + 0.0039*T_airK^2 + 26.352*T_airK + 447.29;
    h_h2o = -6e-17*T_airK^6 + 8e-13*T_airK^5 - 4e-9*T_airK^4 + 8e-6*T_airK^3 - 0.0028*T_airK^2 + 33.367*T_airK + 5.7479;

    % Products enthalpies at T_p
    hp_co2 = -2e-6*T_p^3 + 0.0118*T_p^2 + 34.701*T_p - 2120.1;
    hp_h2o = -6e-17*T_p^6 + 8e-13*T_p^5 - 4e-9*T_p^4 + 8e-6*T_p^3 - 0.0028*T_p^2 + 33.367*T_p + 5.7479;
    hp_n2  = 3e-12*T_p^4 - 5e-7*T_p^3 + 0.0039*T_p^2 + 26.352*T_p + 447.29;
    hp_o2  = -9e-17*T_p^6 + 9e-13*T_p^5 - 3e-9*T_p^4 + 4e-6*T_p^3 + 0.0019*T_p^2 + 27.751*T_p + 126.42;

    % Energy balance equation
    eq = Q_in + W_in + ...
         1 * hf_o + ...
         a * (h_o2 - 8682) + ...
         b * (h_n2 - 8669) + ...
         c * (-241822 + h_h2o - 9904) ...
         == ...
         Q_out + W_out + ...
         d * (-393520 + hp_co2 - 9364) + ...
         e * (-241820 + hp_h2o - 1904) + ...
         f * (hp_n2 - 8669) + ...
         g * (hp_o2 - 8682);

    % Solve equation
    T_p = vpasolve(eq, T_p, [300, 6000]);

end
