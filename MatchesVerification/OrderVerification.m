function [inlLCS]=OrderVerification(Point1,Point2,flag,old_value)
%gets info = 
% collumn_i:Correspondence_i=>[x,y, Geometric information....., order in y,
%order in x
%returns indices of consistent matches
%old = absolutely consistent, new= with thresholds
%old value:pick LCS method
nCorr=size(Point1,2);
if nargin<4    
    old=1;
else
    old=old_value;
end
verbose=true;%for db mesg
if flag
    if old==1
    %Verify X
    [~,index]=sort(Point1(end,:),'ascend');
    SecondSeq=Point2(end,index );
    Indices1=LongestIncreasing(SecondSeq);
    ind_X_Original=index(Indices1);
    %Verify Y
    Okay_x=Point1(end-1,ind_X_Original);
    SecondSeq=Point2(end-1,ind_X_Original);
    [~,index2]=sort(Okay_x,'ascend');
    SecondSeq=SecondSeq(index2);
    Indices2=LongestIncreasing(SecondSeq);
    %Do indices arithmetic to return correspondences of original sequence
    inlLCS=ind_X_Original(index2(Indices2));
    %the next line cancels Yverification and is put here for testing
    %purposes
    %inlLCS=ind_X_Original;
    elseif old==2
        %fixed threshold approach
    [~,index]=sort(Point1(1,:),'ascend');
    SecondSeq=Point2(1,index );
    Indices1=LongestIncreasing(SecondSeq,false);
    ind_X_Original=index(Indices1);
    %Verify Y
    Okay_x=Point1(2,ind_X_Original);
    SecondSeq=Point2(2,ind_X_Original);
    [~,index2]=sort(Okay_x,'ascend');
    SecondSeq=SecondSeq(index2);
    Indices2=LongestIncreasing(SecondSeq,false);
    %Do indices arithmetic to return correspondences of original sequence
    inlLCS=ind_X_Original(index2(Indices2));
    elseif old==3
       %solve and divide
      %calculate matched area based on first image
      %better use a recursive local function
      forced=true;%make sure at least one run

      [~,indexY]=sort(Point1(2,:),'ascend');
      [~,indexX2]=sort(Point2(1,indexY),'ascend');
      Sizey=Point1(2,indexY(end))-Point1(2,indexY(1));   
      verifiedX=callRecursiveLCS(Point1(1,indexY),Point2(1,indexY),...
          indexX2,forced,Sizey,Point1(2,indexY));
      %do sth and pick only verified points
      nX1=Point1(1,indexY(verifiedX));
      nY1=Point1(2,indexY(verifiedX));
      nY2=Point2(2,indexY(verifiedX));
      [~,indexX]=sort(nX1,'ascend'); 
      Sizex=nX1(indexX(end))-nX1(indexX(1));
      [~,indexY2]=sort(nY2(indexX),'ascend');
      verifiedXY=callRecursiveLCS(nY1(indexX),...
          nY2(indexX),indexY2,forced,Sizex,nX1(indexX));
         %combine and return
         inlLCS=indexY(verifiedX(indexX(verifiedXY)));
    else
       %calculate threshold by first-strict run
       [~,index]=sort(Point1(1,:),'ascend');
       SecondSeq=Point2(1,index );
       Indices1=LongestIncreasing(SecondSeq);
       ind_X_Original=index(Indices1);
       %Verify Y
       Okay_x=Point1(2,ind_X_Original);
       SecondSeq=Point2(2,ind_X_Original);
       [~,index2]=sort(Okay_x,'ascend');
       SecondSeq=SecondSeq(index2);
       Indices2=LongestIncreasing(SecondSeq);
       %Do indices arithmetic to return correspondences of original sequence
       FirstLCS=ind_X_Original(index2(Indices2));
       %Evaluate displacement due to misalignment of axes
       keep_orderX=false(1,numel(ind_X_Original));
       keep_orderX(index2(Indices2))=true;
       keep_orderX=ind_X_Original(keep_orderX);
       Erx=arrayfun(@RotXSameZModelFit,...
           Point1(1,keep_orderX(1:end-4)),Point1(1,keep_orderX(2:end-3)),...
           Point1(1,keep_orderX(3:end-2)),Point1(1,keep_orderX(4:end-1)),...
           Point1(1,keep_orderX(5:end)),...
           Point2(1,keep_orderX(1:end-4)),Point2(1,keep_orderX(2:end-3)),...
           Point2(1,keep_orderX(3:end-2)),Point2(1,keep_orderX(4:end-1)),...
           Point2(1,keep_orderX(5:end)));
       Erx=imerode(Erx(~isnan(Erx)),strel([1 1 1]));
       Erx=max(Erx);
%        [MErx,tmp]=max(Erx);
%        Erx(tmp)=-inf;
%        Merx2=max(Erx);
%        if MErx/Merx2>2
%            Erx=Merx2;
%        else
%            Erx=MErx;
%        end
       minY=min(Point2(2,:));
       maxY=max(Point2(2,:));
       maxDY=maxY-minY;
       maxDYv=max(Point2(2,keep_orderX))-min(Point2(2,keep_orderX));
       Erx=Erx*maxDY/maxDYv;
       Ery=arrayfun(@RotXSameZModelFit,...
           Point1(2,FirstLCS(1:end-4)),Point1(2,FirstLCS(2:end-3)),...
           Point1(2,FirstLCS(3:end-2)),Point1(2,FirstLCS(4:end-1)),...
           Point1(2,FirstLCS(5:end)),...
           Point2(2,FirstLCS(1:end-4)),Point2(2,FirstLCS(2:end-3)),...
           Point2(2,FirstLCS(3:end-2)),Point2(2,FirstLCS(4:end-1)),...
           Point2(2,FirstLCS(5:end)));
%        [MErx,tmp]=max(Ery);
%        Ery(tmp)=-inf;
%        Merx2=max(Ery);
%        if MErx/Merx2>2
%            Ery=Merx2;
%        else
%            Ery=MErx;
%        end
       Ery=max(imerode(Ery(~isnan(Ery)),strel([1 1 1])));
       minY=min(Point2(1,:));
       maxY=max(Point2(1,:));
       maxDY=maxY-minY;
       maxDYv=max(Point2(1,keep_orderX))-min(Point2(1,keep_orderX));
       Ery=Ery*maxDY/maxDYv;

  
       if verbose
           fprintf('maxX %.2f maxy %.2f\n',Erx,Ery)
       end
%        warning('Crazy thrs')
       %eval now
       SecondSeq=Point2(1,index );
       Indices1=LongestIncreasing(SecondSeq,false,Erx);
       ind_X_Original=index(Indices1);
       Okay_x=Point1(2,ind_X_Original);
       SecondSeq=Point2(2,ind_X_Original);
       [~,index2]=sort(Okay_x,'ascend');
       SecondSeq=SecondSeq(index2);
       Indices2=LongestIncreasing(SecondSeq,false,Ery);
       inlLCS=ind_X_Original(index2(Indices2));
    end
        
else
    %do nothing
    inlLCS=1:nCorr;
    return
end
end

function VerInd=callRecursiveLCS(Xsort1,X2,index_sort2,forced,size_now,Conj1)
Least_im=200;
T=0.05;%threshold for errors

if (size_now>Least_im) ||(forced)
    %run on those pts
    %LongestIncreasing based on X2 sort, because X1=constrained in location
    %also need to change "old" because it is a true false thingy in LongestINcreasing and it is not so here
    SecondSeq=Xsort1(:,index_sort2);
    Indices=LongestIncreasing(SecondSeq,false,size_now*T);
    if numel(Indices<3)
        VerInd=index_sort2(Indices);
        return;
    end;
    %divide in half and call in each half
    Index_new=index_sort2(Indices);
    Indexoriginal=false(1,size(Xsort1,2));
    Indexoriginal(Index_new)=true;
    Xsort1_new=Xsort1(:,Indexoriginal);
    X2_new=X2(:,Indexoriginal);
    Conj1_new=Conj1(:,Indexoriginal);
    countF=0;
    temp_index=1:size(Xsort1,2);
    for i=1:size(Xsort1,2)
       if~(Indexoriginal(i))
           countF=countF+1;
       else
           temp_index(i)=temp_index(i)-countF;
       end
    end
    index_sort2_new=(temp_index(Index_new));
    size_new=Conj1_new(end)-Conj1_new(1);
    size_new=size_new/2;
    separator=find(Conj1_new>(Conj1_new(1)+size_new),1);
    %call left
    IndicesLeft=callRecursiveLCS(Xsort1_new(1:(separator-1)),...
        X2_new(1:(separator-1)),index_sort2_new(index_sort2_new<separator)...
        ,false,size_new,Conj1_new(1:separator-1));
    %call right
    IndicesRight=callRecursiveLCS(Xsort1_new((separator:end)),...
        X2_new((separator:end)),index_sort2_new(index_sort2_new>=separator)...
       -separator+1 ,false,size_new,Conj1_new(separator:end));
    %combine
    VerInd=1:size(X2,2);
    VerInd=VerInd(Indexoriginal);
    VerInd=[VerInd(IndicesLeft),VerInd(IndicesRight+separator-1)];
    
else
   VerInd=1:size(X2,2); 
end
end

function Er=RotXSameZModelFit(I1,I2,I3,I4,It,II1,II2,II3,II4,IIt)
%I=image in which i sort, II=image in which i find LCS and thus measure
%deviation
%normalise for num stability + equal weight in all elements

Ibar=(I1+I2+I3+I4)/4;
IIbar=(II1+II2+II3+II4)/4;
I1=I1/Ibar-1;
I2=I2/Ibar-1;
I3=I3/Ibar-1;
I4=I4/Ibar-1;
II1=II1/IIbar-1;
II2=II2/IIbar-1;
II3=II3/IIbar-1;
II4=II4/IIbar-1;

A=[ I1 1 -II1*I1 -II1;...
    I2 1 -II2*I2  -II2;...
    I3 1 -II3*I3 -II3;...
    I4 1 -II4*I4 -II4];
[u,s,v]=svd(A);
s(end,end)=0;
Ared=u*s*v';
X=null(Ared);
C=X(3)/(Ibar*IIbar);
D=X(4)/IIbar-C*Ibar;
A=X(1)/Ibar+C*IIbar;
B=X(2)-A*Ibar+IIbar*D+C*Ibar*IIbar;
IIt_est=(A*It+B)/(C*It+D);
Er=abs(IIt_est-IIt);
end



    