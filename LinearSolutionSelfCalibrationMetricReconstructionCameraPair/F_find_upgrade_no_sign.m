function [ H_est1,f1,f2,Pr_1,pinf1,Pr_2,pinf2,H_est2] = F_find_upgrade_no_sign( P_est,F )
%starting from a projective reconstruction, compute focal lengths and p_inf
%INPUT ARGS Pest<- projective matrix, F<-fundamental matrix
%OUTPUT ARGS estimated rectifying homographies, resulting matrices focal
%lengths and pinf coordinates

[S,~,b]=FormSystem(P_est);

%a second round to reverse the image roles
P_est2=vgg_P_from_F(F');
[S2,~,b2]=FormSystem(P_est2);


RR=rref([S,b]);
f1f=sqrt(RR(1,7));
f2f=sqrt(RR(2,7));

%solve second system

RR2=rref([S2,b2]);
f1f2=sqrt(RR2(1,7));
ke_f=f2f/f1f2;
f1=f1f;
f2=f1f2;

%polynomial assuming in third equation, z is the only variable
%two solutions obtained,the right one and another
g1=RR(3,7);
g2=RR(5,7);
g3=RR(4,7);
b=RR(5,6);
g=RR(4,6);
const=g3^2 + g2^2/f1^2-g1;
p2_coef=- 2*f1^2*g*g3- 2*b*g2;
p2_sq_coef= f1^4*g^2+b^2*f1^2+f1^2;
sol=roots([p2_sq_coef,p2_coef,const]);
p2_k=sol(1);
p3_k=g3-g*f1^2*p2_k;
p1_k=(g2-b*f1^2*p2_k)/f1^2;
pinf1=[p1_k,p2_k,p3_k];

%now return the metric reconstruction, the homography restores the
%rotation matrix and gets translation vector up to scale
H_est1=[ke_f^-1*diag([f1,f1,1]),zeros(3,1);-ke_f^-1*[p1_k,p2_k,p3_k]*diag([f1,f1,1]),1];
%H_est1=[diag([f1,f1,1]),zeros(3,1);-[p1_k,p2_k,p3_k]*diag([f1,f1,1]),1];
Pr_1=P_est*H_est1;


    
 

    
    %second solution
    p2_k=sol(2);
    p3_k=g3-g*f1^2*p2_k;
    p1_k=(g2-b*f1^2*p2_k)/f1^2;
    H_est2=[ke_f^-1*diag([f1,f1,1]),zeros(3,1);-ke_f^-1*[p1_k,p2_k,p3_k]*diag([f1,f1,1]),1];
    Pr_2=P_est*H_est2;
    pinf2=[p1_k,p2_k,p3_k];
    

end

