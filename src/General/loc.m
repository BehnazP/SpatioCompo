function [Ndata, index]=loc(location,Data)
% LOC set the matrix of Data same size as location
% 
% [Ndata, index]=loc(location,Data)
% input:
% location      a matrix containg [Lon Lat] that the data should match
% Data          a matrix containing [lon lat Data]       
% 
% output:
% Ndata         a new matrix of data with location [Lon Lat]
% index         an indicator matrix indicating which location of data exist
%               in location
%
% LOC.m 2018-07-15 Behnaz@pirzamanbin.name$

Ndata = nan(size(location,1), size(Data,2));
% set the coordinate of Data be the coordinate of locations
Ndata(:,1:2) = location(:,1:2);

%index matrix 
index = zeros(size(Data,1),1);

% choose values of data if the coordinate is same
for i=1:size(Data,1)
    for j=1:size(location,1)
        if(Data(i,1)==location(j,1) && Data(i,2)==location(j,2)&& size(Data,2)>3)
        Ndata(j,3:end)=Data(i,3:end);
        index(i) = 1;
        elseif(Data(i,1)==location(j,1) && Data(i,2)==location(j,2)&& size(Data,2)==3)
        Ndata(j,3)=Data(i,3);
        index(i) = 1;
        elseif(Data(i,1)==location(j,1) && Data(i,2)==location(j,2)&& size(Data,2)==2)
        index(i) = 1;
        end
    end
end
end
