
function [xn, yn, xf, yf] = generate_trajectories(varargin)
% function generate_trajectories
% The function simulates trajectories of N individuals moving together as a
% rigid body along a randomly generated common trajectory. The positions of
% each individual are not constant with respect to the rigid body, but they
% are affected by time correlated noise
%
% INPUT:
% theta: relative angle of the individuals with respect to the orientation of
% the common trajectory
% nIndividuals: number of individuals in the simulation
% focalRank: rank of the focal individual
% p: other parameters for the trajectories
% spaceUnit: unit of space (e.g. 'metres')
% isSavingFigures: if 1, saves plots of the trajectories as eps graphs
% 
% OUTPUT:
% xn, yn: coordinates of the neighbours (members of the flock other than
% the focal individual)
% xf, yf: coordinates of the focal individual
% 
% Written by:
% Andrea Perna
% http://www.perna.fr
%
% Date:
% 2014 / 04 / 18


if size(varargin, 2) > 0
    theta = varargin{1};
else
    theta = pi/2; % relative angle between the two particles (where the focal particle is with respect to its neighbour)
end

if size(varargin, 2) > 1
    nIndividuals = varargin{2};
else
    nIndividuals = 2; % relative angle between the two particles (where the focal particle is with respect to its neighbour)
end

if size(varargin, 2) > 2
    focalRank = varargin{3};
    if (focalRank > nIndividuals)
        warndlg('the rank of the focal individual cannot be larger than the total number of individuals in the flock', 'error');
        focalRank = randi(nIndividuals);
    end
else
    focalRank = 2; % the focal individual is the first, second or third in a group of three
end

if size(varargin, 2) > 3
    p = varargin{4};
else
    p =[]; % trajectory parameters
end

if size(varargin, 2) > 4
    spaceUnit = varargin{5};
else
    spaceUnit = 'space units'; % relative angle between the two particles (where the focal particle is with respect to its neighbour)
end

if size(varargin, 2) > 5
    isSavingFigures = varargin{6};
else
    isSavingFigures = 0; % relative angle between the two particles (where the focal particle is with respect to its neighbour)
end


if isempty(p)
    % mosquitofish-like : s0 = 5, sV = 0.2, turnR = 0.02, noiseTCorr = 20, Tcutoff = 300
    % pigeon-like : s0 = 5, sV = 0.2, turnR = 0.02, noiseTCorr = 20, or 100, Tcutoff = 300
    
    %% some parameters
    
    % parameters of the trajectory of the neighbour
    p.tL = 2^12; % trajectory length
    p.s0 = 5; % initial speed
    p.sV = 0.2; % speed variation
    p.turnR = 0.02; % max turning rate in radiants per time step
    p.Tcutoff = 300; % this is the period, in time units per cycle of the filter cutoff
    % It means that on average there will be modulations over Tcutoff time
    % units, or smoother
    
    
    % parameters of the mutual position of the focal individual and the
    % neighbour
    p.r = 5; % relative distance between the two particles
    p.noiseTCorr = 20; % temporal correlation of "noise" around the desired position
    % noise can result from sampling noise (e.g. GPS error), or from the
    % incomplete ability of the focal particle to match the target
    p.noiseSpStdT = p.r/2; % spatial standard deviation of the error around the target position (tangential). We can assume that this error is not the same along the two axes tangential and orthogonal to the direction of movement
    p.noiseSpStdO = p.noiseSpStdT; % spatial standard deviation, orthogonal
    
end


% check that speed is not negative
if p.s0 - abs(p.sV) < 0
    errordlg('speed cannot be negative', 'error');
    return
end

%% generate the common trajectory of the flock

rn1 = rand(1, p.tL);
rn2 = rand(1, p.tL);

% remove DC
rn1 = rn1 -mean(rn1);
rn2 = rn2 -mean(rn2);

% fourier transform
r1ft = fft(rn1);
r2ft = fft(rn2);

% create a low-pass filter
sigma=p.tL/p.Tcutoff; % sigma is the standard deviation of the gaussian filter. If we used a step function as a filter instead, the value of sigma would represent the highest frequency of the signal, in cycles per total length (i.e. a value of N means N cycles throughout the whole trajectory)
f = exp(-((((1:p.tL) -p.tL/2 -1).^2)/(2*sigma*sigma))); % create a filter
fShift = ifftshift(f);
fShift(1)=0; % remove the DC (not really necessary because the input signal has no DC)
fShift(end/2) = 0; % this is just for symmetry (because fShift(N) = fShift(end-N+2))

% actually filter the random numbers (now I am using the same filter for both turning and speed; one might say that one of these changes with higher frequency than the other)
r1filt = r1ft .* fShift;
r2filt = r2ft .* fShift;
% reverse Fourier transform (I take the real part to remove approximation errors, but the imaginary part should be really small in this case)
lowPassRN1 = real(ifft(r1filt));
lowPassRN2 = real(ifft(r2filt));

%% Alternative method to generate random trajectories
% lambda = 1;
% [ii,jj] = meshgrid(1:tL);
% K = lambda^2*exp(-(ii-jj).^2/Tcutoff^2);
% C = K; %- K([1, end], :)* inv([K(1,[1, end]); K(end,[1 ,end])]) * K(:, [1, end]);
% lowPassRN1 = mvnrnd(zeros(1,tL), C);
% lowPassRN2 = mvnrnd(zeros(1,tL), C);





%% compute speed and turning angle at each time step
s = p.s0 + lowPassRN1 * p.sV / max(abs(lowPassRN1)); % the speed is equal to the average speed s0 + correlated random oscillations
turn = lowPassRN2 * p.turnR / max(abs(lowPassRN2)); % the turning is given by the random oscillations, scaled so that they do not exceed the maximum turning rate

% plot of speed, acceleration and turning rate
% figure
% subplot(3,1,1)
% plot(s)
% ylabel('speed (U/time step)');
% subplot(3,1,2)
% plot(diff(s))
% ylabel('acc (U/time step)');
% subplot(3,1,3)
% plot(turn)
% ylabel('turn (rad/time step)');

% the direction of movement at time t is the cumulative of all changes of
% direction from 0 to t (assuming the starting direction is 0)
bearing = cumsum(turn);

% x and y coordinates of the neighbour
xn = nan(1, p.tL);
yn = nan(1, p.tL);

% start at the origin
xn(1) = 0;
yn(1) = 0;

% move the neighbour from the starting position, according to speed and bearing
for ii = 1:p.tL-1
    xn(ii+1) = xn(ii) + s(ii)*cos(bearing(ii)); % one might question whether speed and bearing should refer to X(t) - X(t-1) or to X(t+1) - X(t)
    yn(ii+1) = yn(ii) + s(ii)*sin(bearing(ii)); % as these are all random quantities, we don't care too much here.
end





%% generate the noise for the position of particles around the trajectories

% random numbers to generate the noise around the position
grn = randn(2*nIndividuals, p.tL); % gaussian random numbers


% we want these numbers to be correlated in time I repeat a similar Fourier
% filtering operation to the one used above
% fourier transform
grnft = fft(grn, [], 2); % performs the fft transform on each row of length p.tL

% create a low-pass filter
sigma=p.tL/p.noiseTCorr;
f = exp(-((((1:p.tL) -p.tL/2 -1).^2)/(2*sigma*sigma))); % creates a filter
fShift = ifftshift(f);
fShift(1)=0; % remove the DC (not really necessary because the input signal has no DC)
fShift(end/2) = 0; % this is just for symmetry (because fShift(N) = fShift(end-N+2))

% actually filter the random numbers (now I am using the same filter for both turning and speed; one might say that one of these changes with higher frequency than the other)
grnfilt = grnft .* repmat(fShift, 2*nIndividuals, 1);

% reverse Fourier transform (I take the real part to remove approximation errors, but the imaginary part should be really small in this case)
lowPassGNoise = real(ifft(grnfilt,[],2));

% scale standard deviation so that the value matches the desired one
noiseT = lowPassGNoise(1:2:end,:) * p.noiseSpStdT ./ repmat(std(lowPassGNoise(1:2:end,:), 0, 2), 1, p.tL); % tangential noise on the position
noiseO = lowPassGNoise(2:2:end,:) * p.noiseSpStdO ./ repmat(std(lowPassGNoise(2:2:end,:), 0, 2), 1, p.tL); % orthogonal noise on the position


xj = nan(nIndividuals, p.tL);
yj = nan(nIndividuals, p.tL);

individualPositions = linspace(- (nIndividuals-1)/2, (nIndividuals-1)/2, nIndividuals);
for jj = 1:length(individualPositions)
xj(jj, :) = xn + individualPositions(jj) * p.r*cos(bearing + theta) + noiseT(jj,:).*cos(bearing) - noiseO(jj,:).*sin(bearing);
yj(jj,:) = yn + individualPositions(jj) * p.r*sin(bearing + theta) + noiseT(jj,:).*sin(bearing) + noiseO(jj,:).*cos(bearing);
end

xf = xj(focalRank,:); % if rank is 1, the focal individual is at a negative distance relative to bearing + theta (if theta = 0 focal is behind)
yf = yj(focalRank,:);

xn = xj(setdiff(1:nIndividuals, focalRank), :);
yn = yj(setdiff(1:nIndividuals, focalRank), :);

% transform all into column vectors
xn = xn';
yn = yn';
xf = xf';
yf = yf';


colours = 'kcbmgy'; % colors for the trajectories of different neighbours
figure
hold on
for jj = 1:(nIndividuals -1)
plot(xn(:,jj),yn(:,jj), 'Color', colours(mod(jj-1, length(colours))+1));
plot(xn(1:p.Tcutoff:end,jj), yn(1:p.Tcutoff:end,jj), 'Marker', '.', 'LineStyle', 'none', 'Color', colours(mod(jj-1, length(colours))+1));
end
plot(xf,yf, 'r', xf(1:p.Tcutoff:end), yf(1:p.Tcutoff:end), 'r.'); % trajectory of the focal individual in red

axis equal
set(gca, 'TickDir', 'out', 'LineWidth', 1);
set(gca, 'XLim', [min(min(xn))-100, max(max(xn)) + 100], 'YLim', [min(min(yn))-100, max(max(yn)) + 100]);
set(gca, 'FontSize', 16, 'FontName', 'Arial');
xlabel(spaceUnit); ylabel(spaceUnit);
drawnow;
if isSavingFigures
    print(gcf, '-dpng', '-r300', 'trajectory.png');
    print(gcf, '-depsc', '-tiff', '-r300', 'trajectory.eps');
end



