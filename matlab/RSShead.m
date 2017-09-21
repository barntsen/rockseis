function head = RSShead(data)
% RSShead : Creates a header for a RSS formatted file
%
% Call :
% head = RSShead(data);
%
%   Author(s): W. Weibull, 10/09-17
%   Copyright 2017 Wiktor Weibull

MAXDIMS = 9;
head.geometry.N = zeros(MAXDIMS,1);
head.geometry.D = zeros(MAXDIMS,1);
head.geometry.O = zeros(MAXDIMS,1);

Dsz = size(data);

head.data_format = 4;
head.header_format = 4;
head.type = 1;
head.Nheader = 0;
head.Ndims = 1;

for i=1:length(Dsz)
    head.geometry.N(i) = Dsz(i);
end

for i=1:length(Dsz)
    head.geometry.D(i) = 1.0;
end

for i=1:length(Dsz)
    head.geometry.O(i) = 0.0;
end
