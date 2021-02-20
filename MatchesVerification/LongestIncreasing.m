function [inlier_index]=LongestIncreasing(Sequence,old,dx)
%old = strictly increasing
%new=relaxed
%Last mod:Jan28:added dx input
if isempty(Sequence)
    inlier_index=[];
    return
end
if nargin<2
    old=true;
end
nelements=length(Sequence);
top=zeros(nelements+1,1);
index=zeros(nelements,1);
Procesed=zeros(nelements,1);
Rmost=1;
if old
for i=1:nelements
    low=1;
    high=Rmost;
    Picked=Sequence(i);
    Dhl=high-low;
    while (Dhl>1)
       mid= low +floor((high-low)/2);
       check= (Picked<top(mid));
       if check
           high=mid;
       else
           low=mid;
       end
       Dhl=high-low;
    end
    if top(low)>Picked
        place=low;
    else 
        place=high;
    end
    top(place)=Picked;
    index(place)=i;
    if place>1
    Procesed(i)=index(place-1);
    end
    if place==Rmost
        Rmost=Rmost+1;
    end
        
end
else
    if nargin<3
        Dx=120;
    else
    Dx=dx;
    end
for i=1:nelements
    low=1;
    high=Rmost;
    Picked=Sequence(i);
    Dhl=high-low;
    while (Dhl>1)
       mid= low +floor((high-low)/2);
       check= (Picked<top(mid));
       if check
           high=mid;
       else
           low=mid;
       end
       Dhl=high-low;
    end
    if top(low)>Picked
        place=low;
    else 
        place=high;
    end
    top(place)=Picked-Dx;
    index(place)=i;
    if place>1
    Procesed(i)=index(place-1);
    end
    if place==Rmost
        Rmost=Rmost+1;
    end
        
end  
end
    
Rm=Rmost-1;
inlier_index=nan(Rm,1);
last=index(Rm);

for i=Rm:-1:1
    inlier_index(i)=last;
    last=Procesed(last);
    
end

end