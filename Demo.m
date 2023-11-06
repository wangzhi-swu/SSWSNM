% clc
clear
close all

addpath dataset
addpath bootMenu

dataName = 'myDc';
load(dataName);
imOri = myDc;
[m,n,p] = size(imOri);
sigma = [0.050, 0.075, 0.100, 0.125, 0.150, 0.200, 0.300, 0.400];
sig = {'p050','p075','p100','p125','p150','p200','p300', 'p400'};
deltaSets = [0.03,0.03,0.03,0.03,0.03,0.03,0.01,0.01];
lambdaSets = [2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 3.5, 3.5];
p1Sets = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.75, 0.75];

parOneCube.iter = 1;     % 总迭代次数
% parOneCube.delta = 0.03; % 全波段当一块时的最优值 spectral
% parOneCube.lamada = 2.3; % 全波段当一块时的最优值 spectral
parOneCube.c = 2*sqrt(2);     % 全波段当一块时的最优值 spectral
% fprintf('para: delta = %2.2f c = %2.2f ',parOneCube.delta,parOneCube.c); 
psnrBandwise = zeros(length(sigma), p);
ssimBandwise = zeros(length(sigma), p);
fsimBandwise = zeros(length(sigma), p);
errs = zeros(1,50);
results = zeros(3,8);
for noiK = [8:8]
    % for noiK = 1:length(sigma)
    parOneCube.delta = deltaSets(noiK); % 全波段当一块时的最优值 spectral
    parOneCube.lamada = lambdaSets(noiK); % 全波段当一块时的最优值 spectral
    fprintf('para: delta = %2.2f lambda = %2.2f \n',parOneCube.delta,parOneCube.lamada); 
    sigmaMat = sigma(noiK)*ones(1,p);   
    noiName = strcat('case2_', sig{noiK});
    load(noiName, '-mat');
    imNoi = eval(noiName);
   % imNoi = imNoi(1:20,1:20,:);
    psnrIn = myPSNR(imOri, imNoi, 1); 
    ssimIn = mySSIM(imOri, imNoi, 1);
    fprintf('dataName : %s, sigma = %2.4f \n',noiName, sigma(noiK));
    fprintf('input image : psnrIn = %2.4f, ssimIn = %2.4f \n', psnrIn, ssimIn);
    
    estImgAllBand = imNoi;
    estImgAllBand2D = myCubeTo2D(estImgAllBand);
    imNoi2D = myCubeTo2D(imNoi);
    imOri2D = myCubeTo2D(imOri);
    prePSNR = -1;
    norm_x = norm(estImgAllBand2D,'fro');
    TempestImgAllBand2D = estImgAllBand2D; 
    % allBandIter = [4, 6, 10, 14, 17, 28, 32, 56];%case 1
    allBandIter = [4, 7, 11, 13, 21, 28, 32, 56]; %case2
    %% regard all band as a similiar patch and solved by WSNM --- (spectral)
     for one = 1: (allBandIter(noiK))
        estImgAllBand2D = estImgAllBand2D + parOneCube.delta*(imNoi2D - estImgAllBand2D);
        sigmaBand = parOneCube.lamada*sqrt(abs(sigmaMat-mean((imNoi2D - estImgAllBand2D).^2)));
        estImgAllBand2D = onePatEst_wsnm(estImgAllBand2D, sigmaBand, parOneCube, p1Sets(noiK));  
         
        psnrOne = myPSNR(imOri,  reshape(estImgAllBand2D,[m,n,p]), 0);
        ssimOne = mySSIM(imOri,  reshape(estImgAllBand2D,[m,n,p]), 0);
        err = norm(estImgAllBand2D - TempestImgAllBand2D,'fro')^2/((norm_x^2));

         errs(one) = err;
         prePSNR = psnrOne;
         preErr = err;
         % Stop Criteria
%          if(preErr < 1e-4)
%              fprintf('STOP at %d \n',one - 1);
%              break;
%          end
         TempestImgAllBand2D = estImgAllBand2D;
     %   fprintf('iter %d: psnrOne = %2.4f, ssimOne = %2.4f err = %f \n',one, mean(psnrOne), mean(ssimOne), err);
     end
     fprintf('iter %d: Denoise by Step 1: psnrOne = %2.4f, ssimOne = %2.4f err = %f \n',one, mean(psnrOne), mean(ssimOne), err);

    
    
    %% wsnm in spatial domain
    Par = ParSet(imOri);
    Par.nSig0 = sigmaMat;
    Par.I = imOri;
    Par.nim = reshape(estImgAllBand2D,[m,n,p]);
    % Par.Iter 每次迭代的结果
    Par.PSNR  =   zeros( Par.Iter, Par.p);
    Par.SSIM  =   zeros( Par.Iter, Par.p);
    Par.FSIM  =   zeros( Par.Iter, Par.p);
    
    [estImg,Par] = MCWSNM_Denoising(Par.nim, Par.I, Par);
    estImg(estImg>1)=1;
    estImg(estImg<0)=0;
    estImgAllBand = estImg;
    %%
    psnrOut = myPSNR(imOri, estImg, 0);
    ssimOut = mySSIM(imOri, estImg, 0);
    fsimOut = myFSIM(imOri, estImg, 0);
    fprintf('Denoise by Step 2：psnrOut=%2.4f, ssimOut=%2.4f, fsimOut=%2.4f \n', mean(psnrOut),mean(ssimOut),mean(fsimOut));
    % result(noiK) = [mean(psnrOut),mean(ssimOut),mean(fsimOut)]';
   %%
    results(:,noiK) = [mean(psnrOut),mean(ssimOut),mean(fsimOut)]';
    psnrBandwise(noiK,:) = psnrOut;
    ssimBandwise(noiK,:) = ssimOut;
    fsimBandwise(noiK,:) = fsimOut;
    switch noiK
        case 1
            save sswsnmOut2/estImgp050 estImg;
        case 2
            save sswsnmOut2/estImgp075 estImg;
        case 3
            save sswsnmOut2/estImgp100 estImg;
        case 4
            save sswsnmOut2/estImgp125 estImg;
        case 5
            save sswsnmOut2/estImgp150 estImg;
        case 6
            save sswsnmOut2/estImgp200 estImg;
        case 7
            save sswsnmOut2/estImgp300 estImg;
        case 8
            save sswsnmOut2/estImgp400 estImg;
        otherwise
            printf('out of the length of sigma');
    end
    
end

save sswsnmOut2/psnr1-8 psnrBandwise
save sswsnmOut2/ssim1-8 ssimBandwise
save sswsnmOut2/fsim1-8 fsimBandwise
