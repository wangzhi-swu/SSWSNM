function [ Xb ] = solve_wsnm( Y, Cb,K)
iter = 1;
Par.rho = 0.07;
p = 0.9;
Z = zeros(size(Y));
A = zeros(size(Y));
Temp = zeros(size(Y));
X = Y;

while iter < K
    iter = iter + 1;
    Temp1 = X + A;
    [U, SigmaTemp, V] =  svd(full(Temp1), 'econ');
    Temp = diag(SigmaTemp);
    for i=1:1
    W_Vec    =   (Cb)./(Temp.^(1/p)+ eps );
    Z       =   solve_Lp_w(diag(SigmaTemp), W_Vec, p);
    % Temp     =   s1;
    end
    A = Par.rho*(Z - Temp1);
    % Par.rho = min(1e4, Par.mu * Par.rho);
    %     end
end
Xb = Z;
end

