function head = RSShead(data, type)
% RSShead : Creates a header for a RSS formatted file
%
% Call :
% head = RSShead(data, type=1);
%
%   Author(s): W. Weibull, 10/09-17
%   Copyright 2017 Wiktor Weibull

if (nargin < 2)
    head.type=1;
else
    head.type=type;
end

MAXDIMS = 9;
head.geometry.N = zeros(MAXDIMS,1);
head.geometry.D = zeros(MAXDIMS,1);
head.geometry.O = zeros(MAXDIMS,1);

Dsz = size(data);

head.data_format = 4;
head.header_format = 4;

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


if (head.type==2)
    head.Nheader=4;
   for i=1:Dsz(2)
       head.coordinates(i).s(1) = 0;
       head.coordinates(i).s(2) = 0;
       head.coordinates(i).g(1) = 0;
       head.coordinates(i).g(2) = 0;
   end
end

if (head.type==3)
    head.Nheader=6;
   for i=1:Dsz(2)
       head.coordinates(i).s(1) = 0;
       head.coordinates(i).s(2) = 0;
       head.coordinates(i).s(3) = 0;
       head.coordinates(i).g(1) = 0;
       head.coordinates(i).g(2) = 0;
       head.coordinates(i).g(3) = 0;
   end
end
