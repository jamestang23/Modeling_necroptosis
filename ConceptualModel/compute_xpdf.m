%modality testing for marie

function [x2, n, b] = compute_xpdf(x)
  x2 = reshape(x, 1, prod(size(x)));
  %[n, b] = histcounts(x2, 12,'BinLimits',[30,1440+30],'Normalization','pdf');
  [n, b] = histcounts(x2, 24,'BinLimits',[0,24],'Normalization','pdf');
  %[n, b] = histcounts(x2, 12,'BinLimits',[min(min(x)),max(max(x))],'Normalization','pdf');
  %b(1)=[];
  % This is definitely not probability density function
  x2 = sort(x2);
  % downsampling to speed up computations
  %x2 = interp1 (1:length(x2), x2, 1:1000:length(x2));
end

% histogram
% [n, b] = histcounts(x2, 12,'BinLimits',[30,1440+30],'Normalization','pdf');
% bar(n,b);