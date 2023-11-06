function Par = ParSet(imOri)
[Par.m,Par.n,Par.p] = size(imOri);
Par.win       = 40;            % Non-local patch searching window
Par.nlsp      = 70;            % Initial Non-local similar Patch number
Par.Constant  = 2 * sqrt(2);   % Constant num for the weight vector
Par.Innerloop = 2;             % InnerLoop Num of between re-blockmatching
Par.ps        = 6;             % Patch size, larger values would get better performance, but will be slower
Par.step      = 4;             % The step between neighbor patches, smaller values would get better performance, but will be slower
Par.Iter      = 10;            % total iter numbers
Par.display   = true;

Par.maxIter = 10;
Par.model   = 1;
Par.delta   = 0;            % iterative regularization parameter
Par.mu      = 1+eps;
% Par.rho     = 3;              % In final version, this parameter will be changed

% this parameter is not finally determined yet
Par.lambda = 0.8;
% Par.lambda = 8; % set 3 for real
Par.lambda2 = 0;
Par.p2 = 0.9;
