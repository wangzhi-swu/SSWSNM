function  [Z] =  MCWSNM_GST( Y, SigmaRow,SigmaCol, Par )
% Initializing optimization variables
% Intialize the weight matrix W
if Par.lambda2 == 0 
mNSig = min(SigmaRow);
W = (mNSig+eps) ./ (SigmaRow+eps);
else
% W = (mNSig+eps) ./ (SigmaRow+eps);
mNSig = sqrt(min(SigmaRow));
W1 = (mNSig+eps) ./ (sqrt(SigmaRow)+eps);
% W1 = exp( - Par.lambda*mean(SigmaRow, 2));
W2 = 1 ./ (sqrt(SigmaCol) + eps); 
end

% mNSig = min(NSig);
% W = (mNSig+eps) ./ (NSig+eps);
% multi-weighted
if Par.lambda2 == 0 
WY = diag(W) * Y;
tempX = WSNM_2w(WY, Par.Constant,Par.p2,mNSig);
Z = diag(1./W) * tempX;
else
WY = diag(W1) * Y * diag(W2);
tempX = WSNM_2w(WY, Par.Constant,Par.p2,mNSig);
Z = diag(1./W1) * tempX * diag(1./W2);
end
return;
