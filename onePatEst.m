function y = onePatEst(x, sigmaBand, parOneCube)

[U,S,V] =   svd(full(x),'econ'); 
C  = parOneCube.c*sqrt(size(x,2))*2*sigmaBand.^2;

positive = (S+eps).^2-4*repmat(C,size(S,1),1);
ind=find (positive > 0);
valuePositive=length(ind); 

SigmaX = max(S(ind)-eps+sqrt(positive(ind)),0)/2;
y =  U(:,1:valuePositive)*diag(SigmaX)*V(:,1:valuePositive)';
