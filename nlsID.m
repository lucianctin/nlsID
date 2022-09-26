function [tt, yy, resnorm, residual, elim_times, elim_freq] = ...
    nlsID(t,x,tend)

% nlsID is a nonlinear time-domain system identification code. For the full
% description, see Constantin et al (2022) doi.org/10.3390/app12157860
% Copyright (c) 2022 Lucian Constantin, lucian.constantin@bristol.ac.uk

Fs = 1/(t(2)-t(1)); % signal sampling frequency

% get starting frequencies, select them on the FFT
[f, ~] = FindComponents(x, Fs);

Noverlap = NLSprops.Nw - NLSprops.Ns;

% detrend full series (linear)
disp('... linear detrending ...')
% x = detrend(x);
x = x - mean(x(round(end/2):end));


% start ID from maximum
[~,maxind] = max(abs(x));
x = x(maxind:end);
t = t(maxind:end) - t(maxind);

% truncate to tend
x0 = x; t0 = t;
x = x(t<=tend); t = t(t<=tend);

if NLSprops.plotflag
    fig =  figure('Visible', 'off');
    ax1 = axes('Parent', fig);
    hold(ax1, 'on');
    plot(ax1,t0,x0,'k-','LineWidth',2); grid minor
end

% define fit signal
% p(1) - Amplitude
% p(2) - zeta
% p(3) - wn (natural undamped frequency)
% p(4) - phi
p01 = [];

dimprob = length(f);

for i = 1:dimprob
    p01 = [p01 x(1) 0 f(i)*2*pi 0];
end

% define frequency and amplitude bounds
minlim = 1-NLSprops.freqdev; maxlim = 1+NLSprops.freqdev;
minamp = -NLSprops.ampdev; maxamp = NLSprops.ampdev;

% set lower and upper optimization bounds
LB1 = []; UB1 = [];
for i = 1:length(f)
    LB1 = [LB1 minamp*abs(x(1)) -inf f(i)*2*pi*minlim -inf];
    UB1 = [UB1 maxamp*abs(x(1))  inf f(i)*2*pi*maxlim  inf];
end

% set optimization options, tolerances
opts = optimoptions('lsqcurvefit','FiniteDifferenceType','central','display','off','FunctionTolerance',1e-6,'TolFun',1e-6,'TolX',1e-6,'MaxFunctionEvaluations',5e4,'MaxIterations',1e4);

% create structures to store results
tt{1} = []; yy{1} = []; resnorm{1} = []; residual{1} = [];

stopflag1 = 1;
elim_times = []; elim_freq = [];

k1 = 0;
while stopflag1
    % set current window indices
    cr1 = round(k1*(NLSprops.Nw-Noverlap)*1/f(1)*Fs+1:k1*(NLSprops.Nw-Noverlap)*1/f(1)*Fs+1 + NLSprops.Nw*1/f(1)*Fs);
    
    if cr1(end)>length(t)
        break
    end

    k1 = k1+1;

    % extract current signal window
    cursig1 = x(cr1);
    curt1 = t(cr1);
    t10 = curt1(1);

    % time at which signal info is collected is the mean of the current window time
    tt{1}(k1) = t(cr1(1)); 
    curt10 = curt1 - t10;
    % perform least squares fit on current window
    [p1, resnormm, residuall] = lsqcurvefit(@exp_func,p01,curt10(:),cursig1(:),LB1,UB1,opts);

    % verify signal component existence criteria 
    modflag = 0; % flag that signal was modified, means that fit must be done again
    for i = 0:(length(f)-1)
        if p1(i*dimprob+2) > NLSprops.dampthresh
            LB1(i*dimprob+1:i*dimprob+dimprob) = zeros(1,dimprob);
            UB1(i*dimprob+1:i*dimprob+dimprob) = zeros(1,dimprob);
            modflag = 1;
            fprintf('Component %.0f (f = %.1f Hz) eliminated at t = %.2fs\n',i+1,f(i+1),t(cr1(1)));
            elim_times = [elim_times t(cr1(1))]; elim_freq = [elim_freq f(i+1)];
        end
    end
    
    % if any component was eliminated, fit again
    if modflag
        [p1, resnormm, residuall] = lsqcurvefit(@exp_func,p01,curt10(:),cursig1(:),LB1,UB1,opts);
    end
    
    % store results
    for i = 0:(length(f)-1)
        yy{i+1}(k1,:) = p1(i*dimprob+1:i*dimprob+dimprob);
    end

    resnorm{1}(k1) = resnormm;
    residual{1}(k1,:) = residuall;
    
    % set initial conditions for next window
    p01 = p1;

    if NLSprops.plotflag
        sig_reconstruct = exp_func(p1,curt10);
        plot(ax1,curt1,sig_reconstruct,'color',rand(1,3))
        if NLSprops.drawinstant
            set(fig, 'Visible', 'on');
        end

    end

end


