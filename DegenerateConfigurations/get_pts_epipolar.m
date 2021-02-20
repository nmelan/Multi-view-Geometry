function [X_1,X_2]=get_pts_epipolar(P,f1)
%get 100 pts
%uniformly sample epipolar lines, first projection matrix=[k1,0]
X_1=zeros(2,100);
X_2=zeros(2,100);
Pr=P(:,1:3);
v=P(:,4);
for i=1:100
    repeat=1;
    while(repeat)
        try
            repeat=0;
            x_1=-2+4*rand(1);
            y_1=-1+2*rand(1);
            X=diag([1/f1,1/f1,1])*[x_1;y_1;1];
            PsPhOm=Pr*X;
            a_left=(-v(1)-2*v(3))/(PsPhOm(1)+2*PsPhOm(3));
            y_left=(a_left*PsPhOm(2)+v(2))/(a_left*PsPhOm(3)+v(3));
            a_right=(v(2)-v(3))/(PsPhOm(3)-PsPhOm(2));
            x_right=(a_right*PsPhOm(1)+v(1))/(a_right*PsPhOm(3)+v(3));
            lambda=rand(1);
            x_im2=lambda*[-2;y_left]+(1-lambda)*[x_right;1];
        catch
            repeat=1;
        end
    end
        X_1(:,i)=[x_1;y_1];
        X_2(:,i)=[x_im2(1);x_im2(2)];
end