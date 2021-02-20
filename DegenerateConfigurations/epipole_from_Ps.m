% epipole from two projections
function e10 = epipole_from_Ps(P1,P0)
   e10 = zeros(3,1);
   for i=1:3
      e10(i) = det([P1; P0(i,:)]);
   end;
   e10 = e10/norm(e10);
%end;
