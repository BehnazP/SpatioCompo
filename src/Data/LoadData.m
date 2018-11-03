%% load data
load('Data.mat')
load('Data_old.mat')

%choose a time windows 
LPJ = Data_old.LPJ.(time);%dimension of smaller data 
LPJ_Rev = Data_old.LPJ_Rev.(time);%dimension of smaller data 

Elev = Data.Elev;
Elev = Elev(Elev(:,2)>=45,:);
Elev = loc(LPJ(:,1:2),Elev);%dimension of smaller data 

adj = Data.adj.(time);
adj = adj(adj(:,2)>=45,:);
adj = loc(LPJ(:,1:2),adj);%dimension of smaller data 

Elev_Rev = Data.Elev_Rev.(time);
Elev_Rev = Elev_Rev(Elev_Rev(:,2)>=45,:);
Elev_Rev = loc(LPJ_Rev(:,1:2),Elev_Rev);%dimension of smaller data 

adj_Rev = Data.adj_Rev.(time);
adj_Rev = adj_Rev(adj_Rev(:,2)>=45,:);
adj_Rev = loc(LPJ_Rev(:,1:2),adj_Rev);%dimension of smaller data     

Elev0 = Elev;
Elev0(Elev0(:,3)<0,3) = 0;
Elev0_Rev = Elev_Rev;
Elev0_Rev(Elev0_Rev(:,3)<0,3) = 0;

tElev_Rev = (Elev0_Rev(:,3)-mean(Elev0_Rev(:,3)))/std(Elev0_Rev(:,3));
tElev = (Elev0(:,3)-mean(Elev0(:,3)))/std(Elev0(:,3));

Rev = Data.Rev.(time);
Rev = loc(LPJ_Rev(:,1:2),Rev);%dimension of smaller data

[n, m] = size(Rev);
D = m-2;
d = D-1;

[long ,lat] = meshgrid(min(Rev(:,1)):max(Rev(:,1)), ...
              min(Rev(:,2)):max(Rev(:,2)));
sz = size(long);
y_Rev = Rev(:,3:5);

%match REVEALS locations to grid locations and create observation matrix
grid = bsxfun(@minus, long(:), Rev(:,1)').^2 + bsxfun(@minus, lat(:), Rev(:,2)').^2;
[~,I_Rev] = min(grid);
A_REV = sparse(1:length(I_Rev), I_Rev, 1, length(I_Rev), numel(long));

%match LPJ locations to grid locations and create observation matrix
grid = bsxfun(@minus, long(:), LPJ(:,1)').^2 + bsxfun(@minus, lat(:), LPJ(:,2)').^2;
[~,I_LPJ] = min(grid);

A_LPJ = sparse(1:length(I_LPJ), I_LPJ, 1, length(I_LPJ), numel(long));
I_Lpj = reshape(sum(A_LPJ,1)~=0,size(long));

%clean up
clear grid 

%HYDE
load('HYDE_Johan.mat')
HYDE = loc(LPJ(:,1:2),HYDE.LPJ.lu.(time));
if exist('delta','var')==0;delta=1e-5;end
%new LPJ
[LPJnew, ~] = loc(LPJ(:,1:2),Data.LPJ.(time));
adjH = adjust(LPJnew,HYDE);
adjH(:,3:5) = NoZeroOne(adjH(:,3:5),delta);

[adjH_Rev, ~] = loc(Rev(:,1:2),adjH);

%old LPJ
adjH_old = adjust(Data_old.LPJ.(time),HYDE);
adjH_old(:,3:5) = NoZeroOne(adjH_old(:,3:5),delta);

[adjH_old_Rev, ~]=loc(Rev(:,1:2),adjH_old);

if isempty(model)
    error('No model has been selected')
end

switch model
    case 'LPJ_KK10_RCA3'%LPJ-GUESS adjusted with KK10 forced with RCA3
        alradj = Data_old.alradj.(time)(:,3:4);
        alradj_Rev = Data_old.alradj_Rev.(time)(:,3:4);

        p = 4 ;
        B_REV = [ones(length(y_Rev),1) alradj_Rev tElev_Rev];
        B_LPJ = [ones(length(LPJ),1) alradj tElev];
        
    case 'LPJ_KK10_ESM'%LPJ-GUESS adjusted with KK10 forced with ESM
        alradj_Rev = alr(adj_Rev(:,3:5));
        alradj = alr(adj(:,3:5));
        
        p = 4 ;
        B_REV = [ones(length(y_Rev),1) alradj_Rev tElev_Rev];
        B_LPJ = [ones(length(LPJ),1) alradj tElev];
        
    case 'LPJ_HYDE_RCA3'%LPJ-GUESS adjusted with HYDE forced with RCA3
        alradjH = alr(adjH_old(:,3:5));
        alradjH_Rev = alr(adjH_old_Rev(:,3:5));

        p = 4 ;
        B_REV = [ones(length(y_Rev),1) alradjH_Rev tElev_Rev];
        B_LPJ = [ones(length(LPJ),1) alradjH tElev];
        
    case 'LPJ_HYDE_ESM'%LPJ-GUESS adjusted with HYDE forced with ESM
        alradjH_Rev = alr(adjH_Rev(:,3:5));
        alradjH = alr(adjH(:,3:5));

        p = 4 ;
        B_REV = [ones(length(y_Rev),1) alradjH_Rev tElev_Rev];
        B_LPJ = [ones(length(LPJ),1) alradjH tElev];
        
    case 'intercept'
        p = 1 ;
        B_REV = ones(length(y_Rev),1);
        B_LPJ = ones(length(LPJ),1);
        
    case 'elevation'
        p = 2;
        B_REV = [ones(length(y_Rev),1) tElev_Rev];
        B_LPJ = [ones(length(LPJ),1) tElev];
    otherwise
        error('wrong model name')
end

%combined observation/regression matrices
A_x = blkdiag(A_REV,A_REV);
B_Rev = blkdiag(B_REV,B_REV);

