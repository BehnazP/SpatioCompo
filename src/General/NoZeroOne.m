function output = NoZeroOne(data,delta)
% NOZEROONE eleminated zeros and ones from D-compositiona data and instead 
% replace them with delta and 1-delta
% 
% output = NoZeroOne(data,delta)
% input:
% data      data between [0,1]
% delta     small value to be replaced by 0 and 1-delta
%
% output:
% output    new data between (0,1)
%
% NOZEROONE.m 2018-07-15 Behnaz@pirzamanbin.name$

D = size(data,2);

for i = 1:D
    indx0 = find(data(:,i) == 0);
    data(indx0,i) = delta;
    indx1 = find(data(:,i) == 1);
    data(indx1,i) = 1 - delta;
end

if D>1
    output =  bsxfun(@rdivide,data,sum(data,2));
else
    output = data;
end

end