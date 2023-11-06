function  [X] =  WSNM_2w( WY, C, p, NSig)
    
    [U,SigmaWY,V] =   svd(full(WY),'econ');
    PatNum       = size(WY,2);
    s = diag(SigmaWY);
    Temp         =   sqrt(max( s.^2 - PatNum*NSig^2, 0 ));
    % Temp = s;
    TempC = (C*sqrt(PatNum)*NSig^2);
    s1 = zeros(size(s));
  
    for i=1:4
        W_Vec    =   TempC./( Temp.^(1/p) + eps );
        % W = W_Vec ./ diag(SigmaW);
        % Weight vector
        %W_Vec    =   (C*sqrt(PatNum)*NSig^2)./( Temp + eps );
       	s1       =   solve_Lp_w(s, W_Vec, p);
       	Temp     =   s1;
    end
    SigmaX = diag(s1);
    X =  U*SigmaX*V'; 
return;
