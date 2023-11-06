function p = myPSNR(img1, img2, type)
% ==========================================
% function p = myPSNR(img1, img2, type)
% type = 0, 返回各个波段的psnr
% type = 1, 返回总的平均psnr
% ==========================================

if max(img1(:)) > 10
    img1 = img1 / max(img1(:));
end
if max(img2(:)) > 10
    img2 = img2 / max(img1(:));
end

[h, w, ch] = size(img1);
img1 = reshape(img1, [h*w ch]);
img2 = reshape(img2, [h*w ch]);
mse = mean( (img1-img2).^2 );
psnr = 10 * log10( 1 ./ mse );
if type == 0
    p = psnr;
elseif type ==1
    p = mean( psnr );
end