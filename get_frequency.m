function get_frequency(varargin)
% Computes FFT
x = varargin{1};
Fs = varargin{2};
fig = varargin{3};

L = length(x);
f1 = Fs*(0:(L/2))/L; % frequency values
Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
figure(fig); loglog(f1,P1,'-');
grid minor;
end