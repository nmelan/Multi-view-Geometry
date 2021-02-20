function [ R ] = rotate_align( v1,v2 )
%R*v1=v2
axis=vgg_contreps(v1);
axis=axis*v2;
axis=axis/norm(axis);
if ~isfinite(axis(1))
    axis=vgg_contreps(v1)*rand(3,1);
    axis=axis/norm(axis);
end
v1n=v1/norm(v1);
v2n=v2/norm(v2);
cos_th=v1n'*v2n;
th=acos(cos_th);
R=rotation_angle_axis(th,axis);


end

