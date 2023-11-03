function [maxDirCorrTime, maxDirCorrValue, maxSpeedCorrTime, maxSpeedCorrValue, secondPeakFocalTime, secondPeakFocalValue, secondPeakNeighbourTime, secondPeakNeighbourValue] = speed_and_directional_correlation(sf, sn, bearingf, bearingn, varargin)
% function speed_and_directional_correlation
% The function computes correlations between the speed and direction of a
% focal individual and the speed and direction of its neighbour(s), based
% on Pearson's correlation coefficient.
%
% INPUT:
% sf: speed sequence of the focal individual
% sn: speed sequence of the neighbour(s)
% bearingf: bearing sequence of the focal individual
% bearingn: bearing sequence of the neighbour(s)
% 
% Written by:
% Andrea Perna
% http://www.perna.fr
%
% Date:
% 2014 / 04 / 18


if size(varargin, 2) > 0
    maxLag = varargin{1};
else
    maxLag = 100; % maximum lag for speed and directional correlation, in time steps
end


%% Speed correlations

% temporal autocorrelations of speed
[rF, lagsF] = xcorr(sf - mean(sf), sf - mean(sf), maxLag, 'coeff');
[valMaxF,iMaxF,valMinF,iMinF] = extrema(rF); % extrema of correlation

% figure, plot(lagsF,rF)
% hold on;
% plot(lagsF(iMaxF),valMaxF,'g.'); % local maxima
% plot(lagsF(iMinF),valMinF,'r.') % local minima

if length(iMaxF) >=2
    secondPeakFocalTime = lagsF(iMaxF(2));
    secondPeakFocalValue = valMaxF(2);
else
    secondPeakFocalTime = NaN;
    secondPeakFocalValue = NaN;
end

for jj = 1:size(sn,2)
[rN(:,jj), lagsN(:,jj)] = xcorr(sn(:,jj) - mean(sn(:,jj)), sn(:,jj) - mean(sn(:,jj)), maxLag, 'coeff');
[valMaxN,iMaxN,valMinN,iMinN] = extrema(rN(:,jj)); % extrema of correlation

% figure, plot(lagsN(:,jj),rN(:,jj))
% hold on;
% plot(lagsN(iMaxN),valMaxN,'g.'); % local maxima
% plot(lagsN(iMinN),valMinN,'r.') % local minima

if length(iMaxN) >=2
    secondPeakNeighbourTime = lagsN(iMaxN(2));
    secondPeakNeighbourValue = valMaxN(2);
else
    secondPeakNeighbourTime = NaN;
    secondPeakNeighbourValue = NaN;
end


% % temporal correlation of speed between focal individual and neighbour
[r(:,jj), lags(:,jj)] = xcorr(sn(:,jj) - mean(sn(:,jj)), sf - mean(sf), maxLag, 'coeff');
% If there are missing values, consider using something like the lines
% below instead of xcorr(...)
% counter = 0;
% for tau = -maxLag:1:maxLag
%     counter = counter + 1;
%     lags(counter) = tau;
%     if tau < 0
%         [r(counter), pval(counter)] = corr(sf(-tau+1:end)', sn(1:end+tau)');
%     else
%         [r(counter), pval(counter)] = corr(sf(1:end-tau)', sn(tau +1:end)');
%     end
% end

% figure, plot(lags,r)
[maxSpeedCorrValue(jj), I] = max(r(:,jj));
maxSpeedCorrTime(jj) = lags(I,jj);

%% Directional correlations

clear lagsDir rDir;
counter = 0;
for tau = -maxLag:1:maxLag
    counter = counter + 1;
    lagsDir(counter, jj) = tau;
    if tau < 0
        rDir(counter, jj) = mean(cos(bearingf(-tau+1:end) - bearingn(1:end+tau, jj)));
    else
        rDir(counter, jj) = mean(cos(bearingf(1:end-tau) - bearingn(tau +1:end, jj)));
    end
end


% figure, plot(lagsDir,rDir)
[maxDirCorrValue(jj), I] = max(rDir(:,jj));
maxDirCorrTime(jj) = lagsDir(I,jj);


end

