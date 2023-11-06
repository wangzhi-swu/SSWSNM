function y = myAddNoise(x, par, type)
% ==============================================
% function y = myAddNoise(x, par, type)
% x:      original data
% par:    added noise level
%          -- sigma: type=0 (��������׼���������,��������)
%          -- snr:   type=1 (����ֵ����ȼ�����)
% type:   noise added type
% y:      output(noisy) data
% ==============================================

[m,n,p] = size(x);
y=zeros(m,n,p);

if max(x(:))>10
    x = x/max(x(:));
end

if type == 0
    for i = 1:p
        y(:,:,i) = x(:,:,i) + par(i)*randn(m,n); % ����imnoise
    end
elseif type ==1
    for i=1:p
        band=x(:,:,i);
        varNoise= norm(band(:))^2/(length(band(:)) * (10^ (par/10)));
        y(:,:,i)=band+sqrt(varNoise)*randn(size(band));
    end
end


