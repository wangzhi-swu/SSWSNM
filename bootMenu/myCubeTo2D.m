function y = myCubeTo2D (x)
[m,n,p] = size(x);
y = zeros(m*n,p);
for i = 1:p 
    temp = x(:,:,i);
    y(:,i) = temp(:);
end