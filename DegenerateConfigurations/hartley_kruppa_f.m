% [fs,cond] = hartley_kruppa_f(F,e,pp0,pp1)
% Given a fundamental matrix and one epipole derived from cameras
% with zero skew, equal aspect ratio, and known principal points 
% pp0,pp1, find the unknown focal lengths f0,f1 using Hartley's
% Kruppa equation based method.

function [fs,cond] = hartley_kruppa_f(F10,e10,pp0,pp1)
   if (nargin<4)
      if (nargin<3) pp0=[0;0;1]; end; 
      pp1=pp0; 
   end;
   me = mcross(e10);
   Fp = F10' * pp1;
   mep = me * pp0;
   FIF = F10' * diag([1,1,0]) * F10;
   FpF = Fp * Fp';
   eIe = me * diag([1,1,0]) * me';
   epe = mep * mep';
   A1 = zeros(6,2);
   A2 = A1;
   for i=1:3
      for j = 1:i
	 A1(i*(i-1)/2+j,:) = [FIF(i,j), FpF(i,j)];
	 A2(i*(i-1)/2+j,:) = [eIe(i,j), epe(i,j)];
      end;
   end;
   % [F10,e10],pp1'*F10*pp0
   % scaling didn't improve results:
   %   A1 = A1/norm(A1,1);
   %   A2 = A2/norm(A2,1);
   % [A1,A2]
   [U,S,V] = svd([A1,A2]);
%diag(S)'
   sol = V(:,4);
   cond = [S(4,4)/S(3,3),S(3,3)/S(1,1)];
   fs = sqrt([sol(3)/sol(4), sol(1)/sol(2)]);
%end;
