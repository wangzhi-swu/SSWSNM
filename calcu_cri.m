%% quality assess
addpath(genpath('quality_assess'));
[~,~,p] = size(myDc);
npsnr_v = zeros(1,p);nssim_v = zeros(1,p);nfism_v = zeros(1,p);nsam_v = zeros(1,p);nergas_v = zeros(1,p);
mpsnr_v = zeros(1,p);mssim_v = zeros(1,p);mfism_v = zeros(1,p);msam_v = zeros(1,p);mergas_v = zeros(1,p);
for i=1:1:p
    J=255*myDc(:,:,i);
    K=255*case1_p300(:,:,i);
    I=255*estImg(:,:,i);
    [npsnr_v(i),nssim_v(i),nfsim_v(i),nergas_v(i),nmsam_v(i)] = MSIQA(J,K);
    [mpsnr_v(i),mssim_v(i),mfsim_v(i),mergas_v(i),mmsam_v(i)] = MSIQA(J,I);
end
fprintf('psnr = %.4f , ssim = %.4f, fsim = %.4f , ergas = %.4f \n',mean(mpsnr_v),mean(mssim_v),mean(mfsim_v),mean(mergas_v));



