%% Independent MATLAB reproduction and pulse-bin sensitivity study
% Units: rate = microgray/hour, time = seconds, dose = microgray.
% All interval rates are exact wall-time means. No instrument inference.
% MATLAB R2019b+ (script with local functions).
clear; clc; close all;
outdir = fileparts(mfilename('fullpath'));
if isempty(outdir), outdir = pwd; end
TOL = 1e-11;

%% 1. Independently reproduce published synthetic controls
starts = [0:100:800, 900:20:980];
ends = [100:100:900, 920:20:1000];
assert(numel(starts)==14 && numel(ends)==14);
constant = 2*ones(size(starts));
adaptive = [ones(1,9), 10*ones(1,5)];
linearEdges = [0 100 250 500 1000];
linearRates = 1 + 0.002*(linearEdges(1:end-1)+linearEdges(2:end))/2;
controlNames = {'constant';'linear';'adaptive_pulse';'narrow_pulse'};
truth = [5/9; 5/9; 19/36; 5/18];
weighted = [sum(constant.*(ends-starts))/3600; ...
            sum(linearRates.*diff(linearEdges))/3600; ...
            sum(adaptive.*(ends-starts))/3600; ...
            sum([10 10].*[50 50])/3600];
unweighted = [mean(constant)*1000/3600; ...
              mean(linearRates)*1000/3600; ...
              mean(adaptive)*1000/3600; NaN];
assert(all(abs(weighted-truth)<TOL));
assert(abs(unweighted(2)-0.4652777777777778)<TOL);
assert(abs(unweighted(3)-1.1706349206349207)<TOL);
controlTable = table(controlNames,truth,weighted,unweighted, ...
    'VariableNames',{'Case','AnalyticDose_uGy','WeightedDose_uGy','UnweightedDose_uGy'});
writetable(controlTable,fullfile(outdir,'control_reproduction.csv'));
disp(controlTable);

%% 2. Bin width x pulse onset x window end = 5 x 3 x 3 = 45 cases
widths = [5 10 20 50 100];
onsets = [40 45 50];
windowEnds = [48 52 100];
rows = zeros(numel(widths)*numel(onsets)*numel(windowEnds),11);
k=0;
for p=onsets
    for w=widths
        edges = 0:w:100;
        assert(edges(end)==100);
        binLeft=edges(1:end-1); binRight=edges(2:end);
        % Exact integrated dose in each bin; pulse occupies [p,p+10).
        pulseOverlap = max(0,min(binRight,p+10)-max(binLeft,p));
        meanRate = 100*pulseOverlap./(binRight-binLeft);
        fullDose = sum(meanRate.*(binRight-binLeft))/3600;
        assert(abs(fullDose-5/18)<TOL);
        for tEnd=windowEnds
            % True dose on [0,tEnd): overlap with the actual pulse.
            exactDose = 100*max(0,min(tEnd,p+10)-max(0,p))/3600;
            % Proportional allocation assumes constant rate inside cut bin.
            binWindowOverlap = max(0,min(binRight,tEnd)-max(binLeft,0));
            estimate = sum(meanRate.*binWindowOverlap)/3600;
            % Nonnegative-only sharp bounds from full-bin means.
            complete = binWindowOverlap >= (binRight-binLeft)-TOL;
            intersected = binWindowOverlap > 0;
            lower = sum(meanRate(complete).*(binRight(complete)-binLeft(complete)))/3600;
            upper = sum(meanRate(intersected).*(binRight(intersected)-binLeft(intersected)))/3600;
            assert(exactDose >= lower-TOL && exactDose <= upper+TOL);
            err = estimate-exactDose;
            if exactDose > 0
                pct = 100*err/exactDose;
            else
                pct = NaN; % Relative error undefined for zero ground truth.
            end
            k=k+1;
            rows(k,:) = [p,w,tEnd,fullDose,exactDose,estimate,err,abs(err),pct,lower,upper];
        end
    end
end
assert(k==45);
results = array2table(rows,'VariableNames',{'PulseStart_s','BinWidth_s','WindowEnd_s', ...
    'FullDose_uGy','TrueWindowDose_uGy','EstimatedWindowDose_uGy', ...
    'SignedError_uGy','AbsoluteError_uGy','RelativeError_pct', ...
    'NonnegativeLowerBound_uGy','NonnegativeUpperBound_uGy'});
writetable(results,fullfile(outdir,'pulse_bin_error_table.csv'));
disp(results);

%% 3. Plot absolute error versus bin width (partial windows only)
fig = figure('Color','w','Position',[100 100 950 560]);
styles = {'-o','-s','-^'};
for idx=1:numel(onsets)
    p=onsets(idx);
    subplot(1,2,1); hold on;
    mask = results.PulseStart_s==p & results.WindowEnd_s==48;
    plot(results.BinWidth_s(mask),results.AbsoluteError_uGy(mask),styles{idx}, ...
        'LineWidth',1.6,'MarkerSize',6,'DisplayName',sprintf('Start %d s',p));
    subplot(1,2,2); hold on;
    mask = results.PulseStart_s==p & results.WindowEnd_s==52;
    plot(results.BinWidth_s(mask),results.AbsoluteError_uGy(mask),styles{idx}, ...
        'LineWidth',1.6,'MarkerSize',6,'DisplayName',sprintf('Start %d s',p));
end
for idx=1:2
    subplot(1,2,idx); grid on;
    xlabel('Uniform bin width (s)'); ylabel('Absolute dose error (\muGy)');
    windowLabels = [48 52];
    title(sprintf('Window [0,%d) s',windowLabels(idx)));
    legend('Location','best');
end
sgtitle('Pulse timing and bin-width sensitivity (proportional allocation)');
exportgraphics(fig,fullfile(outdir,'pulse_bin_error.png'),'Resolution',200);

%% 4. Additional invariants
assert(all(abs(results.FullDose_uGy-5/18)<TOL));
fullMask=results.WindowEnd_s==100;
assert(all(abs(results.SignedError_uGy(fullMask))<TOL));
% 45-second pulse, 50-second bins, 48-second cutoff: 2/15 estimate, 1/12 truth.
mask=results.PulseStart_s==45 & results.BinWidth_s==50 & results.WindowEnd_s==48;
assert(nnz(mask)==1 && abs(results.EstimatedWindowDose_uGy(mask)-2/15)<TOL);
assert(abs(results.TrueWindowDose_uGy(mask)-1/12)<TOL);
fprintf('PASS: four controls, 45 scenarios, all conservation/boundary checks.\n');
