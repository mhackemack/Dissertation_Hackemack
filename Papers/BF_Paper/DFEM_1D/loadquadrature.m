function loadquadrature(sn)

global snq

% get gauss-legendre abscissae and weights
[a,b]=GLNodeWt(sn);
% put data in snq struct
snq.sn=sn;
snq.mu=a';
snq.w =b';

% sum of the weights
snq.sw = sum(snq.w);

return
end