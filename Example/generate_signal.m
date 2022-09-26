function x = generate_signal(t,p)
% generate test function using p form
x = zeros(size(t));
for i = 1:size(p,2)/4
    N = (i-1)*4;
    x = x + (p(:,N+1)'.*exp(-p(:,N+2)'.*p(:,N+3)'.*t).*cos(p(:,N+3)'.*sqrt(1-p(:,N+2)'.^2).*t+p(:,N+4)'));
end
end


