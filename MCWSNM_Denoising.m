function [rI, Par]   =  MCWSNM_Denoising( nI, I, Par )
rI           = nI;         % Estimated Image
[h, w, ch]  = size(rI);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );
% noisy image to patch
NoiPat =	Image2Patch( nI, Par );  %ÿ�����ε�С�鰴������һ�У�ÿ�д�Сps^2*3
Par.TolN = size(NoiPat, 2); % ÿ��ͨ�����е�С��ĸ���������ÿ��ͨ��ת��һ���ˣ������������С�������ÿ��ͨ����С�������ͬ
Sigma_arrCh = zeros(Par.ch, Par.TolN); % ����ͨ����ÿ��ͨ����С���sigma��ɵľ���
Par.Sigma = sqrt(mean(Par.nSig0.^2));
prePSNR = -1;
tempRI = nI;
for iter = 1 : Par.Iter
    Par.iter = iter;
    % iterative regularization
    rI =	rI + Par.delta * (nI - rI);
    % image to patch
    CurPat =	Image2Patch( rI, Par ); % ȫ������һ�飬����ͨ�����һ��
    % CurPat�ߴ磺 (ps*ps*3) * [(r-ps+1)*(c-ps+1)], Ҳ����˵ÿ����������ͨ����С�����ɵģ��ܹ�С�������
    % ��NoiPat��С��ͬ
    % estimate local noise variance
    for c = 1:Par.ch
        if(iter == 1)
            TempSigma_arrCh = sqrt(   max(  0, repmat(Par.nSig0(c)^2, 1, size(CurPat, 2))... 
                            - mean( (NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) ... % ����ͼ�͵�ǰͼ��ÿ��ͨ����С��������������������ƽ����Ȼ��ȡ��ֵ
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






