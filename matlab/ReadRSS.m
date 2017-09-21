function [data, head] = ReadRSS(filename)
% ReadRSS : Reads a RSS formatted file
%
% Call :
% [Data,Header]=ReadRSS(filename);
% Example:
% [Data, Header] = ReadRSS('../Test/MPI/Acoustic2D/Pshot.rss');
% ot = Header.geometry.O(1); % Time origin
% dt = Header.geometry.D(1); % Sampling interval 
% nt = Header.geometry.N(1); % Number of time samples
% ntr = Header.geometry.N(2); % Number of traces
% SrcCoords = [Header.coordinates(1:ntr).s]; % Source coordinates (x,y,z)
% RecCoords = [Header.coordinates(1:ntr).g]; % Receiver coordinates (x,y,z)
% % Display data:
% t = (0:nt-1)*dt + ot; 
% x = RecCoords(1:Header.Nheader/2:end);
% figure, imagesc(1:ntr, t, Data); 
% set(gca, 'xticklabel', x(get(gca,'xtick')));
%
%   Author(s): W. Weibull, 10/09-17
%   Copyright 2017 Wiktor Weibull

MAGICNUMBER='R0CKS'; %rss identifier
MAGICNUMBERLENGTH=5;
MAXDIMS=9;
data=[];
head=[];

fileID = fopen(filename, 'r');
% Check if file has rss identifier
buff = fread(fileID,MAGICNUMBERLENGTH,'*char')';
if(strcmp(buff, MAGICNUMBER) == false)
    display('Not a valid rss file!');
    return;
end

% Reading header
head.data_format = fread(fileID, 1,'*int');
head.header_format = fread(fileID, 1,'*int');
head.type = fread(fileID, 1,'*int');
head.Nheader = fread(fileID, 1,'*uint64');
head.Ndims = fread(fileID, 1,'*uint64');

for i=1:MAXDIMS
    head.geometry.N(i) = fread(fileID, 1,'*uint64');
end

for i=1:MAXDIMS
    head.geometry.D(i) = fread(fileID, 1,'*float64');
end

for i=1:MAXDIMS
    head.geometry.O(i) = fread(fileID, 1,'*float64');
end

% Reading traces and coordinates
if(head.Nheader)
    data=zeros(head.geometry.N(1), head.geometry.N(2));
    for i=1:head.geometry.N(2)
        for j=1:head.Nheader/2
            head.coordinates(i).s(j) = fread(fileID, 1,'*float32');
        end
        for j=1:head.Nheader/2
            head.coordinates(i).g(j) = fread(fileID, 1,'*float32');
        end
        
        data(:,i) = fread(fileID, head.geometry.N(1),'*float32');
    end
else
    for i=1:MAXDIMS
        if (head.geometry.N(i) > 0)
            Dsz(i) = head.geometry.N(i);
        end
    end
    data = fread(fileID, prod(Dsz),'*float32');
    data = reshape(data, Dsz);
end

fclose(fileID);

end
