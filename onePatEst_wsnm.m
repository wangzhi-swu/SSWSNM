function y = onePatEst_wsnm(x, sigmaBand, parOneCube,p1)
[U,S,V] =   svd(full(x),'econ'); 
C  = parOneCube.c*sqrt(size(x,2))*2*sigmaBand.^2;
% C = parOneCube.c*sqrt(size(x,2))*2;
Temp = diag(S);
% Temp = sqrt(max( diag(S).^2 - size(x,2)*sigmaBand'.^2, 0 ));
s = diag(S);
% p = 0.8;
for i=1:1
    W_Vec    =   (C')./(Temp.^(1/p1)+ eps );
    s1       =   solve_Lp_w1(s, W_Vec, p1);
    Temp     =   s1;
end
SigmaX = diag(s1);
y =  U*SigmaX*V' ;

% positive = (S+eps).^2-4*repmat(C,size(S,1),1);
% ind=find (positive > 0);
% valuePositive=length(ind); 
% SigmaX = max(S(ind)-eps+sqrt(positive(ind)),0)/2;
% y =  U(:,1:valuePositive)*diag(SigmaX)*V(:,1:valuePositive)';
