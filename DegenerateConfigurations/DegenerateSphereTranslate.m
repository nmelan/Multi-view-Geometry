function [f1_erMe,f2_erMe,f1_erKr,f2_erKr]=DegenerateSphereTranslate(mag,NL)
%test near degenerate


%rot axis parallel or perp to translation, then fix angle to pi/2

f1_gt=0.500+1.200*rand(1);
f2_gt=0.500+1.200*rand(1);
Sph=UniformSphereSampling(2);
Pos1=Sph(:,1)*mag;
Pos2=Sph(:,2);
dir1=-Sph(:,1);
dir2=-Sph(:,2);
R1=rotate_align(dir1,[0;0;1]);
R2=rotate_align(dir2,[0;0;1]);
K2=diag([f2_gt,f2_gt,1]);
P2=K2*[R2,-R2*Pos2];
P2=P2*[R1'*diag([1/f1_gt,1/f1_gt,1]),Pos1;zeros(1,3),1];
%solution
[X_1,X_2]=get_pts_epipolar(P2,f1_gt);
X_1=X_1+NL*X_1.*randn(2,100);
X_2=X_2+NL*X_2.*randn(2,100);
F = estimateFundamentalMatrix(X_1',X_2');
P_est=vgg_P_from_F(F);
[~,rnk,~]=FormSystem(P_est);
if (rnk==5)
[ ~,f1,f2,~,~,~,~,~]=F_find_upgrade_no_sign(P_est,F);
f1_erMe=df_eval_error(f1,f1_gt);
f2_erMe=df_eval_error(f2,f2_gt);
else
f1_erMe=NaN;
f2_erMe=NaN;
end
runbougnoux;%output Fb,Fh
if (imag(Fb(1))~=0 || imag(Fb(2))~=0)
    f1_erKr=NaN;
    f2_erKr=NaN;
else
   f1_erKr= df_eval_error(Fb(1),f1_gt);
   f2_erKr= df_eval_error(Fb(2),f2_gt);
end
    