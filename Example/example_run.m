clear

% generate test signal with 4 frequency components and added gaussian noise
rng('default');

Fs = 5000;
tend = 8;
t = linspace(0,tend,Fs*tend);
f = [2 5 20 50]; % [Hz]
zeta = [4 4 5 8]/100;
A = [10 5 3 8];
phi = [0 rand(1,3)*2*pi];

p = [];
for i = 1:length(f)
    p = [p A(i) zeta(i) f(i)*2*pi phi(i)];
end
xref = generate_signal(t,p);

% corrupt signal with white Gaussian noise
x = awgn(xref,30,'measured');

t1 = t; x1 = -sign(x(1))*x;

% plot test signal
figure; plot(t1,x1,'k-');

%%
% apply the nlsID on test signal and see how it performs
[tt, pres1, resnorm, residual] = ...
    nlsID(t1,x1,tend);

%%
% plot results
tnorm = tt{1}*f(1);

figure
for i = 1:length(f)
    plot(tnorm,abs(pres1{i}(:,1)),'o-'); hold on;
end
title('Amplitudes')
xlabel('t/T1'); ylabel('Component amplitude [-]'); set(gca,'FontSize',18);

figure
for i = 1:length(f)
    plot(tnorm,pres1{i}(:,2)*100,'o-'); hold on;
end
xlabel('t/T1'); ylabel('Component damping [%]');  set(gca,'FontSize',18);
title('Damping')

figure
for i = 1:length(f)
    plot(abs(pres1{i}(:,1)),pres1{i}(:,2)*100,'o-'); hold on;
end
xlabel('Component amplitude [-]'); ylabel('Component damping [%]');  set(gca,'FontSize',18);
title('Damping')


figure
for i = 1:length(f)
    plot(tnorm,pres1{i}(:,3)/2/pi,'o-'); hold on;
end
title('Frequencies')
xlabel('t/T1'); ylabel('Component frequency [Hz]');  set(gca,'FontSize',18);

figure
for i = 1:length(f)
    plot(abs(pres1{i}(:,1)),pres1{i}(:,3)/2/pi,'o-'); hold on;
end
title('Frequencies')
xlabel('Component amplitude [-]'); ylabel('Component frequency [Hz]');  set(gca,'FontSize',18);


figure
plot(tnorm,resnorm{1},'o-');
title('2-norm of residuals')
xlabel('t/T1'); ylabel('resnorm [-]');  set(gca,'FontSize',18);


figure
plot(tnorm,residual{1}(:,1:3),'o-');
title('residuals')
xlabel('t/T1'); ylabel('residual [-]'); set(gca,'FontSize',18);


figure 
for i = 1:length(f)
    plot(tnorm,mod(pres1{i}(:,4),2*pi),'o-'); hold on;
end
title('Phase angles')
xlabel('t/T1'); ylabel('phi modulo 2\pi [rad]');  set(gca,'FontSize',18);






