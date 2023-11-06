
function  [X] =  WSNM( Y,  mNSig,  Par)
%     p = Par.pAd;
p = 0.9;
[U,SigmaY,V] =   svd(full(Y),'econ');
[PatNum ]      = size(Y,2);
Temp = diag(SigmaY) ;
% Temp =  sqrt(max( diag(SigmaY).^2 - PatNum*mSig^2, 0 ));
s = diag(SigmaY);
s1 = zeros(size(s));
TempC = Par.Constant * sqrt(PatNum) * mNSig^2;

for i=1:4 
    W_Vec    =   (TempC)./(Temp.^(1/p)+ eps );
    s1       =   solve_Lp_w(s, W_Vec, p);
    Temp     =   s1;
end

%     positiveNum=find (s1 > 0);
%   ind=length(positiveNum);

SigmaX = diag(s1);
X =  U*SigmaX*V' ;
return;
