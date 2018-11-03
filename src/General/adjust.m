function adj = adjust(data1,data2)
% ADJ is a function that adjust the value of data1 (3-composition) with 
% a vector of data2
% 
% adj = adjust(data1,data2)
% input:
% data1     a 3-compositional data to be adjusted with data2; 
%           [Lon, Lat, 3-composition]
% data2     a vector for adjustment of data1; [Lon, Lat, data2]
%
% output:
% adj       new data1 adjusted with the vector data2
%
% NOZEROONE.m 2018-07-15 Behnaz@pirzamanbin.name$

adj(:,1:2) = data1(:,1:2);

for i = 1:3
    if i ~= 3
        adj(:,i+2) = (1-data2(:,3)).*data1(:,i+2);
    else
        adj(:,i+2) = (1-data2(:,3)).*data1(:,i+2) + data2(:,3);
    end
end
    
end