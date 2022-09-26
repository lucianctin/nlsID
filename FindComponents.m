function [f, tlims] = FindComponents(y, Fs)
% Find linear approximation of frequencies
fig = figure;
get_frequency(y,Fs,fig);
hold on;
disp('Select frequency components and hit enter')
[f,tlims] = ginput;
end