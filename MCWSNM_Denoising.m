function [rI, Par]   =  MCWSNM_Denoising( nI, I, Par )
rI           = nI;         % Estimated Image
[h, w, ch]  = size(rI);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );
% noisy image to patch
NoiPat =	Image2Patch( nI, Par );  %每个波段的小块按列拉成一列，每列大小ps^2*3
Par.TolN = size(NoiPat, 2); % 每个通道所有的小块的个数，现在每个通道转成一列了，所以这里的总小块个数和每个通道的小块个数相同
Sigma_arrCh = zeros(Par.ch, Par.TolN); % 三个通道，每个通道的小块的sigma组成的矩阵
Par.Sigma = sqrt(mean(Par.nSig0.^2));
prePSNR = -1;
tempRI = nI;
for iter = 1 : Par.Iter
    Par.iter = iter;
    % iterative regularization
    rI =	rI + Par.delta * (nI - rI);
    % image to patch
    CurPat =	Image2Patch( rI, Par ); % 全部遍历一遍，三个通道组成一列
    % CurPat尺寸： (ps*ps*3) * [(r-ps+1)*(c-ps+1)], 也就是说每个列是三个通道的小块拉成的，总共小块个数列
    % 和NoiPat大小相同
    % estimate local noise variance
    for c = 1:Par.ch
        if(iter == 1)
            TempSigma_arrCh = sqrt(   max(  0, repmat(Par.nSig0(c)^2, 1, size(CurPat, 2))... 
                            - mean( (NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) ... % 噪声图和当前图像每个通道的小块提出来，做差，各个像素平方，然后取均值
                                   - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2 )  )   );  
            %             TempSigma_arrCh = sqrt(abs(repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
        else
            TempSigma_arrCh = Par.lambda*sqrt(   max(  0, repmat(Par.nSig0(c)^2, 1, size(CurPat, 2))...
                            - mean( (NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :)... 
                                   - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2 )  )   );
            %             TempSigma_arrCh = Par.lambda*sqrt(abs(repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
        end
        Sigma_arrCh((c-1)*Par.ps2+1:c*Par.ps2, :) = repmat(TempSigma_arrCh, [Par.ps2, 1]);
    end
    
    %Estimated Local Noise Level
    Sigma_Col = Par.lambda2 * sqrt(abs(repmat(Par.Sigma^2, 1, size(CurPat,2)) - mean((NoiPat - CurPat).^2))); 

    if (mod(iter-1, Par.Innerloop) == 0)
        Par.nlsp = Par.nlsp - 10;  % Lower Noise level, less NL patches
        NL_mat  =  Block_Matching(CurPat, Par); % Caculate Non-local similar patches for each
    end
    % Denoising by MCWSNM
    [Y_hat, W_hat]  =  MCWSNM_Estimation( NL_mat, Sigma_arrCh,Sigma_Col,CurPat, Par );
    rI = PGs2Image(Y_hat, W_hat, Par);
    rI(rI>1)=1;
    rI(rI<0)=0;
    PSNR  = myPSNR( I, rI, 0 );
    SSIM  = mySSIM( I, rI, 0 );
    FSIM  = myFSIM( I, rI, 0 );
    if prePSNR > mean(PSNR)
        rI = tempRI;
        PSNR  = myPSNR( I, rI, 0 );
        SSIM  = mySSIM( I, rI, 0 );
        FSIM  = myFSIM( I, rI, 0 );
   %     fprintf( 'Final Iter = %2.0f, PSNR = %2.4f, SSIM = %2.4f, FSIM = %2.4f \n', iter-1, mean(PSNR), mean(SSIM), mean(FSIM) );
        Par.PSNR(iter-1,:)  =   PSNR;
        Par.SSIM(iter-1,:)  =   SSIM;
        Par.FSIM(iter-1,:)  =   FSIM;
        break;
    end
    prePSNR = mean(PSNR);
    tempRI = rI;    
   % fprintf( 'Iter = %2.0f, PSNR = %2.4f, SSIM = %2.4f, FSIM = %2.4f \n', iter, mean(PSNR), mean(SSIM), mean(FSIM) );
    Par.PSNR(iter,:)  =   PSNR;
    Par.SSIM(iter,:)  =   SSIM;
    Par.FSIM(iter,:)  =   FSIM;
end






