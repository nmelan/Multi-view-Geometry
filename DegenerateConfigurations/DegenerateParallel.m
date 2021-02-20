function [f1_erMe,f2_erMe,f1_erKr,f2_erKr]=DegenerateParallel(thetain,NL)
%test near degenerate


%rot axis parallel or perp to translation, then fix angle to pi/2
axis_rot=rand(1,3);
axis_rot=axis_rot./norm(axis_rot);
theta=thetain;
R=rotation_angle_axis(theta,axis_rot);
f1_gt=0.500+1.200*rand(1);
f2_gt=0.500+1.200*rand(1);
P1=[diag([f1_gt,f1_gt,1]),zeros(3,1)];
%parallel

%90deg case
%end 90 deg case
%
tr=axis_rot;
P2=diag([f2_gt,f2_gt,1])*[R,tr'];
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
    