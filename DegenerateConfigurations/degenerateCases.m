%test degeneracies
%uncomment proper section to test
%


%rot axis parallel or perp to translation, then fix angle to pi/2
axis_rot=rand(1,3);
axis_rot=axis_rot./norm(axis_rot);
%theta=rand(1)*pi/2;
theta=pi/2;
R=rotation_angle_axis(theta,axis_rot);
f1_gt=0.500+1.200*rand(1);
f2_gt=0.500+1.200*rand(1);
P1=[diag([f1_gt,f1_gt,1]),zeros(3,1)];

%parallel

%90deg case
%Kruppa fails, i have rank deficiency(check what can be done)
%end 90 deg case
%
% tr=axis_rot;
% P2=diag([f2_gt,f2_gt,1])*[R,tr'];

%perp
%okay for both

% tr=vgg_contreps(axis_rot)*rand(3,1);
% tr=tr./norm(tr);
% P2=diag([f2_gt,f2_gt,1])*[R,tr];

%on sphere pointing at center
%Kruppa fails, i loose rank

% Sph=UniformSphereSampling(2);
% dir1=-Sph(:,1);
% dir2=-Sph(:,2);
% R1=rotate_align(dir1,[0;0;1]);
% R2=rotate_align(dir2,[0;0;1]);
% K1=diag([f1_gt,f1_gt,1]);
% K2=diag([f2_gt,f2_gt,1]);
% Pos1=-dir1;
% Pos2=-dir2;
% P1=K1*[R1,-R1*Pos1];
% P2=K2*[R2,-R2*Pos2];



%bmvc96
%use xy-plane as common plane

%case1
%kruppa fail, i drop rank
% t=rand(1);
% Vectors=[rand(2,4);t,t,zeros(1,2)];
% R1=rotate_align(Vectors(:,3),[0;0;1]);
% R2=rotate_align(Vectors(:,4),[0;0;1]);
% K1=diag([f1_gt,f1_gt,1]);
% K2=diag([f2_gt,f2_gt,1]);
% Pos1=Vectors(:,1);
% Pos2=Vectors(:,2);
% P1=K1*[R1,-R1*Pos1];
% P2=K2*[R2,-R2*Pos2];

%case2
%kruppa fails i loose rank
t=rand(1);
Vectors=[rand(2,3);t,t,zeros(1,1)];
R1=rotate_align(Vectors(:,3),[0;0;1]);
base=Vectors(:,1)-Vectors(:,2);
ratio=base(2)/base(1);
optical21=rand(1);
optical22=optical21*ratio;
optical23=rand(1);
optical2=[optical21;optical22;optical23];
R2=rotate_align(optical2,[0;0;1]);
K1=diag([f1_gt,f1_gt,1]);
K2=diag([f2_gt,f2_gt,1]);
Pos1=Vectors(:,1);
Pos2=Vectors(:,2);
P1=K1*[R1,-R1*Pos1];
P2=K2*[R2,-R2*Pos2];


%solution
F=vgg_F_from_P(P1,P2);
P_est=vgg_P_from_F(F);
[ ~,f1,f2,~,~,~,~,~]=F_find_upgrade_no_sign(P_est,F);
f1_erMe=df_eval_error(f1,f1_gt);
f2_erMe=df_eval_error(f2,f2_gt);
runbougnoux;