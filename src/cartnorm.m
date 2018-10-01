function[x] = cartnorm(A)
%convenience function - probably need to optimise for numeric stability
[m,n] = size(A);

x = zeros(m,1);
for ii = 1:n
    x = x + A(:,ii).^2;
end
x = sqrt(x);