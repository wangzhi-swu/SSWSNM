addpath(genpath('sswsnmOut2'));
addpath(genpath('dataIn2'));
sig = {'p050','p075','p100','p125','p150','p200','p300', 'p400'};
load('sswsnmOut\myDc.mat');
for noiK = 1:8
denoiName = strcat('sswsnmOut2\estImg', sig{noiK});
load(denoiName, '-mat');
noiName = strcat('dataIn2\case2_', sig{noiK});
load(noiName,'-mat');
output_image = estImg;
% 获取MAT文件中包含的所有变量的名称
var_names = who('-file', noiName);
% 如果MAT文件中只包含一个变量，则加载该变量并将其分配给MATLAB变量
if numel(var_names) == 1
    data = load(noiName,'-mat');
    % 给变量指定新名称
    new_var_name = 'Noisy_Img';
    assignin('base', new_var_name, data.(var_names{1}));
end
%% quality assess
addpath(genpath('quality_assess'));
[~,~,p] = size(myDc);
npsnr_v = zeros(1,p);nssim_v = zeros(1,p);nfism_v = zeros(1,p);nsam_v = zeros(1,p);nergas_v = zeros(1,p);
mpsnr_v = zeros(1,p);mssim_v = zeros(1,p);mfism_v = zeros(1,p);msam_v = zeros(1,p);mergas_v = zeros(1,p);
for i=1:1:p
    J=255*myDc(:,:,i);
    K=255*Noisy_Img(:,:,i);
    I=255*output_image(:,:,i);
    [npsnr_v(i),nssim_v(i),nfsim_v(i),nergas_v(i),nmsam_v(i)] = MSIQA(J,K);
    [mpsnr_v(i),mssim_v(i),mfsim_v(i),mergas_v(i),mmsam_v(i)] = MSIQA(J,I);
end
fprintf('psnr = %.4f , ssim = %.4f, fsim = %.4f , ergas = %.4f \n',mean(npsnr_v),mean(nssim_v),mean(nfsim_v),mean(nergas_v));
fprintf('psnr = %.4f , ssim = %.4f, fsim = %.4f , ergas = %.4f \n',mean(mpsnr_v),mean(mssim_v),mean(mfsim_v),mean(mergas_v));
results(:,noiK) = [mean(mpsnr_v),mean(mssim_v),mean(mfsim_v),mean(mergas_v)]';
end