function res = exp_func(p,tt)
% damped cosines from p structure
res = zeros(size(tt));
for i = 1:length(p)/4
    N = (i-1)*4;
    res = res + (p(N+1)*exp(-p(N+2)*p(N+3)*tt).*cos(p(N+3)*sqrt(1-p(N+2)^2)*tt+p(N+4)));
end
end

