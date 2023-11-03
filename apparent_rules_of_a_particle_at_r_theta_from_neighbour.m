
function [] = apparent_rules_of_a_particle_at_r_theta_from_neighbour(varargin)
% function apparent_rules_of_a_particle_at_r_theta_from_neighbour
% this function generates 120 random trajectories of two (or more)
% particles moving together at a fixed position relative to each other
%
% INPUT:
% attractorDirection
% This is the direction of the target position for an individual relative
% to its neighbour. attractorDirection = 0 produces "flocks" where all
% individuals move on the same segment, tangent to the main trajectory,
% one after the other; attractorDirection = pi/2 produces "flocks" where
% all individuals move side by side on the same segment, orthogonal to the
% main trajectory
% noiseTCorr
% this gives the temporal scale for correlation in the noise around the
% position of each individual. Larger values indicate that the noise
% remains correlated for longer time (an individual observed to the left or
% to the right of its "real" position is likely to be observed on the same
% side for a longer time)
% nIndividuals
% number of individuals in the simulation
%
% Written by:
% Andrea Perna
% http://www.perna.fr
%
% Date:
% 2014 / 04 / 18






if size(varargin, 2) > 0
    attractorDirection = varargin{1};
else
    attractorDirection = pi/2; % relative angle between the two particles (where the focal particle is with respect to its neighbour)
    %0 or pi/2% this expresses the direction of the attractor,
    % that is the position where the focal individual aims at being relative to
    % its neighbour. The position of the focal individual at time t is based on
    % the orientation of the neighbour between time t and time t+1 (they start
    % the path between two time steps being aligned.
    % When the value of theta = 0, the target position is randomly
    % chosen in front or behind the neighbour; when the value is pi/2, the
    % target position of the focal individual is on the side of its neighbour.
end

if size(varargin, 2) > 1
    trajectoryParameters.noiseTCorr = varargin{2};
else
    trajectoryParameters.noiseTCorr = 20;
end

if size(varargin, 2) > 2
    nIndividuals = varargin{3};
else
    nIndividuals = 2;
end


% generate multiple trajectories and store all data together in the
% variable allData

% parameters
isSavingTrajectories = 0;
spaceUnitsForTrajectories = 'm';

% mosquitofish : s0 = 5, sV = 0.2, turnR = 0.02, noiseTCorr = 20, Tcutoff = 300 theta = 0;
% pigeon : s0 = 5, sV = 0.2, turnR = 0.02, noiseTCorr = 20, or 100, Tcutoff = 300 theta = pi/2

% parameters of the trajectory of the neighbour
trajectoryParameters.tL = 2^12; % trajectory length
trajectoryParameters.s0 = 5; % initial speed
trajectoryParameters.sV = 0.2; % speed variation
trajectoryParameters.turnR = 0.02; % max turning rate in radiants per time step
trajectoryParameters.Tcutoff = 300; % this is the period, in time units per cycle of the filter cutoff
% It means that on average there will be modulations over Tcutoff time
% units, or smoother


% parameters of the mutual position of the focal individual and the
% neighbour
trajectoryParameters.r = 5; % relative distance between the two particles
% this is an input parameter: trajectoryParameters.noiseTCorr = 20; % temporal correlation of "noise" around the desired position
% noise can result from sampling noise (e.g. GPS error), or from the
% incomplete ability of the focal particle to match the target
trajectoryParameters.noiseSpStdT = trajectoryParameters.r/2; % spatial standard deviation of the error around the target position (tangential). We can assume that this error is not the same along the two axes tangential and orthogonal to the direction of movement
trajectoryParameters.noiseSpStdO = trajectoryParameters.noiseSpStdT; % spatial standard deviation, orthogonal

allConfigurationClasses = [];
allMaxDirCorrTimes = [];
allMaxDirCorrValues = [];
allMaxSpeedCorrTimes = [];
allMaxSpeedCorrValues = [];
allSecondPeakFocalTimes = [];
allSecondPeakFocalValues = [];
allSecondPeakNeighbourTimes = [];
allSecondPeakNeighbourValues = [];

allFocalData =[];
allNeighbourData = [];
for jj = 1:120
    jj
    close all;
    % generate trajectories
    focalRank = randi(nIndividuals);
    allConfigurationClasses = [allConfigurationClasses; focalRank];
    [xn, yn, xf, yf] = generate_trajectories(attractorDirection, nIndividuals, focalRank, trajectoryParameters, spaceUnitsForTrajectories, isSavingTrajectories);
    
    isSavingTrajectories = 0; % save only the first trajectory
    
    
    % compute bearing and speed as a function of time
    [bearingn, sn] = cart2pol(diff(xn,1,1), diff(yn,1,1));
    [bearingf, sf] = cart2pol(diff(xf,1,1), diff(yf,1,1));
    
    
    % compute acceleration and turning at each time step
    an = diff(sn,1,1);
    af = diff(sf);
    turnn = mod (pi + diff(bearingn,1,1), 2*pi) - pi;
    turnf = mod (pi + diff(bearingf,1,1), 2*pi) - pi;
    
    
    % compute r (distance), theta (relative direction) and phi (relative orientation) of the neighbour
    [theta, r] = cart2pol(xn - repmat(xf, 1, size(xn,2)), yn - repmat(yf, 1, size(xn,2))); % these repmats are useless here because there is only one neighbour (and it would not work with more neighbours below)
    theta = mod (pi + theta(2:end,:) - repmat(bearingf, 1, size(xn,2)), 2*pi) - pi;
    phi = mod(pi + bearingn - repmat(bearingf, 1, size(xn,2)), 2*pi) - pi;
    
    
    % clip the data so that xn(t), sn(t) and an(t) represent respectively
    % the position at time t, the speed between time t-1 and time t, and the
    % acceleration between the speed at time (t-1 to t) and the speed at time
    % (t to t+1)
    
    xn = xn(2:end-1,:);
    yn = yn(2:end-1,:);
    xf = xf(2:end-1);
    yf = yf(2:end-1);
    r = r(2:end-1,:);
    
    bearingn = bearingn(1:end-1,:);
    bearingf = bearingf(1:end-1);
    sn = sn(1:end-1,:);
    sf = sf(1:end-1);
    
    theta = theta(1:end-1,:);
    phi = phi(1:end-1,:);
    
    
    allFocalData = [allFocalData; ...
        xf ...      (1)
        yf ...      (2)
        sf ...      (3)
        af ...      (4)
        turnf ...   (5)
        ];
    
    allNeighbourData = [allNeighbourData; ...
        xn ...      (nNeighbours*0 + 1: nNeighbours*1)
        yn ...      (nNeighbours*1 + 1: nNeighbours*2)
        sn ...      (nNeighbours*2 + 1: nNeighbours*3)
        r ...       (nNeighbours*3 + 1: nNeighbours*4)
        theta ...   (nNeighbours*4 + 1: nNeighbours*5)
        phi ...     (nNeighbours*5 + 1: nNeighbours*6)
        an ...      (nNeighbours*6 + 1: nNeighbours*7)
        turnn...    (nNeighbours*7 + 1: nNeighbours*8)
        ];
    
    % compute directional and speed correlations
    [maxDirCorrTime, maxDirCorrValue, maxSpeedCorrTime, maxSpeedCorrValue, ...
        secondPeakFocalTime, secondPeakFocalValue, secondPeakNeighbourTime, secondPeakNeighbourValue] = ...
        speed_and_directional_correlation(sf, sn, bearingf, bearingn);
    
    % record speed and directional correlation for all trajectories
    allMaxDirCorrTimes = [allMaxDirCorrTimes; maxDirCorrTime];
    allMaxDirCorrValues = [allMaxDirCorrValues; maxDirCorrValue];
    allMaxSpeedCorrTimes = [allMaxSpeedCorrTimes; maxSpeedCorrTime];
    allMaxSpeedCorrValues = [allMaxSpeedCorrValues; maxSpeedCorrValue];
    allSecondPeakFocalTimes = [allSecondPeakFocalTimes; secondPeakFocalTime];
    allSecondPeakFocalValues = [allSecondPeakFocalValues; secondPeakFocalValue];
    allSecondPeakNeighbourTimes = [allSecondPeakNeighbourTimes; secondPeakNeighbourTime];
    allSecondPeakNeighbourValues = [allSecondPeakNeighbourValues; secondPeakNeighbourValue];
    
    
end

isSavingDirCorrDelayFigures = 0;
% plot directional correlation and speed correlation results
% positive values mean the focal individual moved later
plot_ranks_vs_classes(allMaxDirCorrTimes', allConfigurationClasses', 1, isSavingDirCorrDelayFigures, 'trajectory number', 'Anticipation / delay (time steps)', 'focal position', 'dir. corr. delay (positive values: focal is leading)')


% clear all variables except allData
clearvars -except allFocalData allNeighbourData nIndividuals


sf = allFocalData(:,3);
af = allFocalData(:,4);
turnf = allFocalData(:,5);

nNeighbours = nIndividuals - 1;
sn = allNeighbourData(:,nNeighbours*2 + 1: nNeighbours*3);
r = allNeighbourData(:,nNeighbours*3 + 1: nNeighbours*4);
theta = allNeighbourData(:,nNeighbours*4 + 1: nNeighbours*5);
phi = allNeighbourData(:,nNeighbours*5 + 1: nNeighbours*6);


isSavingFigures = 0;
useAllPhiBins = 0; % if 0 it only plots some phi bins (corresponding to nearly aligned particles)
framesPerTimeUnit = 5; % all the parameters are in units per frame. Here I tell how many frames are in a time unit
timeUnit = 's';
spaceUnit = 'm';


make_colourful_plots(sf, sn, r, theta, phi, af, turnf, framesPerTimeUnit, timeUnit, spaceUnit, useAllPhiBins, isSavingFigures);













