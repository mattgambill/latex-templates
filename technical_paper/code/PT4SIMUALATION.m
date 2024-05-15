%{
ECE 5013 Radar Project
Pt4 Simulation Runs the radar
Author: Matt Gambill
Date: 04/21/20
%}
clear, clc;
fprintf('Script Started at %s\n', datestr(now,'HH:MM:SS.FFF'))
tic;
%% Parameters for Target
x0 = 25;            % R0 is initially 30 meters
theta_intial = -10; % theta = azimuth angle, initially at -10 degrees
theta_final = 10;   % theta = azimuth angle, finally at 10 degrees
v0 = [0 10];             % vertical  velocity is 10 m/s

initialTargetPos = [x0 x0*tand(theta_intial)];
finalTargetPos = [x0 x0*tand(theta_final)];
simulationTime = (finalTargetPos - initialTargetPos)/v0;

Np = 64; % Pulses per CPI
fp = 1e3; % pulse sample rate


timePerCPI = Np/fp;
numCPIs = ceil(simulationTime/timePerCPI)+1;

t = linspace(0,simulationTime,numCPIs);
s = [v0(1)*t+initialTargetPos(1); v0(2)*t+initialTargetPos(2)];
actualRange = sqrt(s(1,:).^2+s(2,:).^2);
actualAz = atan(s(2,:)./s(1,:));
%polarplot(actualAz,actualRange)

%%
simRange = zeros(1,numCPIs);
simVel = zeros(1,numCPIs);
simAz = zeros(1, numCPIs);
%% Processing
parfor CPI = 1:numCPIs
    
    [simRange(CPI), simVel(CPI), simAz(CPI)] = radarSimulator(actualRange(CPI),actualAz(CPI),v0,Np,fp,CPI);
    
end
%% Post Processing

velocityError = (abs(simVel - v0(2))/v0(2))';
rangeError = (abs(simRange-actualRange)./actualRange)';
azimuthError = (abs(simAz - actualAz)./actualAz)';
tableData = [rangeError azimuthError velocityError];
saveNormalizedErrorTable

scriptTime = toc/60;

fprintf('The script took %f minutes\n',scriptTime)