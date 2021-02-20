[jcc, cc] = function Jccfun(imnamesN, F_estimates)
%input arguments: number of images, f. length 
Confident_cell=cell(imnamesN,1);
JConfident_cell=cell(imnamesN,1);
for i=1:imnamesN
    for j=1:imnamesN
        if all(isnan(F_estimates{i,j}))
            F_estimates{i,j}=[];
        end
    end
end
F_length=zeros(imnamesN,imnamesN);
for i=1:imnamesN
    for j=1:imnamesN
        F_length(i,j)=length(F_estimates{i,j});
    end
end


Num_points=max(Num_points,Num_points');
N_weight=ones(imnamesN,imnamesN);

ff11=cell(imnamesN,1);
for i=1:imnamesN
    ff11{i}.f1=cell2mat(F_estimates(i,:));
    ff11{i}.f2=cell2mat(F_estimates(:,i)');
    n2=length(ff11{i}.f1);
    ff11{i}.source=nan(n2,1);
    start=1;
    for j=1:imnamesN
        add=F_length(i,j);
        ff11{i}.source(start:start+add-1)=j*ones(add,1);
        start=start+add;
    end
    
    confident_count=zeros(1,n2);
    [ff11{i}.f1,index]=sort(ff11{i}.f1);
    ff11{i}.f2=ff11{i}.f2(index);
    ff11{i}.source=ff11{i}.source(index);
    threshold=0.1;
    low=1;high=1;nh=2;
    for j=1:n2
        picked=ff11{i}.f1(j);
        while (low<=n2)&&(df_eval_error(ff11{i}.f1(low),picked)>threshold)
            low=low+1;
        end
        while(high<n2)&&(df_eval_error(ff11{i}.f1(nh),picked)<threshold)
            high=high+1;
            nh=nh+1;
        end
        %confident_count(j)=nh-low;
        %this is now the weighted sum, way i implement sum BAD
        confident_count(j)=sum(N_weight(i,ff11{i}.source(low:high)));
    end
    Confident_cell{i}=confident_count/n2;
end
for i=1:imnamesN
    confident=Confident_cell{i};
    flocal=ff11{i};
    n2=length(confident);
    match_set=find(F_length(i,:)>0);
    for j=1:n2
        f2now=flocal.f2(j);
        im2=flocal.source(j);
        [~,ind_joint]=ismember(f2now,ff11{im2}.f1);
        conf_now=Confident_cell{im2}(ind_joint);
        confident(j)=confident(j)+conf_now;
        f1now=flocal.f1(j);
        matchset=match_set(match_set~=flocal.source(j));
        for l=matchset
            f1_insecond=ff11{l}.f2;
            er_f1_insecond= arrayfun(@(x) df_eval_error(x,f1now),f1_insecond);
            f_keep=Confident_cell{l}(er_f1_insecond<0.1);
            if ~isempty(f_keep)
                confident(j)=confident(j)+nanmedian(f_keep);
            end
        end
    end
    plot(flocal.f1,confident);
    f_est_joint=flocal.f1(confident>=0.95*max(confident));
    k_cf=median(f_est_joint);
	JConfident_cell{i}=confident/n2;
end


jcc = JConfident_cell;
cc =  Confident_cell;
end