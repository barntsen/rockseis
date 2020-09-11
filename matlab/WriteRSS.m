function WriteRSS(filename, data, head)
% WriteRSS : Writes a RSS formatted file
%
% Call :
% WriteRSS(filename, data, head);
% The head can be created using head = RSShead(data);
% And the modifying the header parameters, such as sampling intervals and coordinates 
% at will. 
%
%   Author(s): W. Weibull, 10/09-17
%   Copyright 2017 Wiktor Weibull

if(nargin < 3) 
    display('Not enough input arguments.');
    return;
end
MAGICNUMBER='R0CKS'; %rss identifier
MAXDIMS=9;

fileID = fopen(filename, 'w');
% Check if file has rss identifier
fwrite(fileID, MAGICNUMBER,'*char');


% Writting header
fwrite(fileID, head.data_format,'*int');
fwrite(fileID, head.header_format,'*int');
fwrite(fileID, head.type,'*int');
fwrite(fileID, head.Nheader,'*uint64');
fwrite(fileID, head.Ndims,'*uint64');

for i=1:MAXDIMS
    fwrite(fileID, head.geometry.N(i),'*uint64');
end

for i=1:MAXDIMS
    fwrite(fileID, head.geometry.D(i), '*float64');
end

for i=1:MAXDIMS
    fwrite(fileID, head.geometry.O(i),'*float64');
end
    
% Writing traces and coordinates
if(head.Nheader)    
    for i=1:head.geometry.N(2)
        for j=1:head.Nheader/2
            fwrite(fileID, head.coordinates(i).s(j),'float32');
        end
        for j=1:head.Nheader/2
            fwrite(fileID, head.coordinates(i).g(j),'float32');
        end
        
        fwrite(fileID, data(:,i),'float32');
    end
else    
    fwrite(fileID, data,'*float32');    
end

fclose(fileID);

end
