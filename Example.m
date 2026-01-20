clc; clear; close all;
%%
% --- Usage Example Script ---

% 1. Initialize all mole fractions to 0
myMoleFractions = zeros(1, 17);

% 2. Define specific mole fractions
myMoleFractions(1)  = 0.5;   % Water
myMoleFractions(2)  = 0.5;   % Methanol
myMoleFractions(3)  = 0;     % Ethylene Glycol
myMoleFractions(4)  = 0;     % 1,2-Propanediol
myMoleFractions(5)  = 0;     % 1,3-Butanediol
myMoleFractions(6)  = 0;     % Diethylene Glycol
myMoleFractions(7)  = 0;     % N,N-Dimethylethanolamine
myMoleFractions(8)  = 0;     % 2-Methyl-2,4-Pentanediol
myMoleFractions(9)  = 0;     % Triethylene glycol (TER)
myMoleFractions(10) = 0;     % Diethylene Glycol Monobutyl Ether
myMoleFractions(11) = 0;     % Glycerol
myMoleFractions(12) = 0;     % DMSO
myMoleFractions(13) = 0;     % 3-Methoxy-1,2-Propanediol
myMoleFractions(14) = 0;     % 1,2,4 butanetriol
myMoleFractions(15) = 0;     % 1,3-Propanediol
myMoleFractions(16) = 0;     % 2,3-Butanediol
myMoleFractions(17) = 0;     % 1,4-Butanediol

% 3. Define VFT setting (true/false) and output filename
useVFT = true;
outputFile = 'MyResults.xlsx';

%Temperature Range
    LowT = -20;
    HighT = 35;
    StepT = 0.5;
%% 
% Run prediction with VFT enabled, saving to 'MyResults.xlsx'
finalResults = ViscosityPredictor(myMoleFractions, LowT, StepT, HighT, useVFT, outputFile);
% Temperature = finalResults(:,1);
% Viscosity=finalResults(:,2);

Temperature = finalResults.Temperature; 
Viscosity   = finalResults.Viscosity;
%% plot
figure;
plot(Temperature, Viscosity, 'b-', 'LineWidth', 2);
grid on;
xlabel('Temperature (°C)');
ylabel('Viscosity (cP)');
title('Predicted Viscosity vs Temperature');
