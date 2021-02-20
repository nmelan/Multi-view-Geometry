% fs = bougnoux_kruppa_f(F10,e10,e01,pp0,pp1)
% Given a fundamental matrix and one epipole derived from cameras
% with zero skew, equal aspect ratio, and known principal points 
% pp0,pp1, find the unknown focal lengths f0,f1 using Bougnoux'
% Kruppa equation based method.

function fs = bougnoux_kruppa_f(F10,e10,e01,pp0,pp1)
   if (nargin<5) 
      if (nargin<4) 
          pp0=[0;0;1];
      end; 
      pp1=pp0; 
   end;
   me0 = mcross(e10);
   me1 = mcross(e01);
   d110 = diag([1,1,0]);
   pFp = pp1' * F10 * pp0;
   pmf0 = pp1'*me1*d110*F10;
   pmf1 = pp0'*me0*d110*F10';
   fs = sqrt(-pFp*[ ...
	  (pmf0*pp0)/(pmf0*d110*F10'*pp1), ...
	  (pmf1*pp1)/(pmf1*d110*F10*pp0)]);
end
