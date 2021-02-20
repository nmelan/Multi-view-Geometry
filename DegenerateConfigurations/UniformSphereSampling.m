function [samples]=UniformSphereSampling(n)
%returns n samples (3xn matrix) on a sphere with radius one
phi=2*pi*rand(n,1);
cosine_theta=-1+2*rand(n,1);
temp=[cos(phi),sin(phi),ones(n,1)];
sin_theta=1-cosine_theta.*cosine_theta;
sin_theta=sqrt(sin_theta);
temp_theta=[sin_theta,sin_theta,cosine_theta];
sampless=temp_theta.*temp;
samples=sampless';
end
