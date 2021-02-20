F10=F;
P_est1=[eye(3),zeros(3,1)];
e10 = epipole_from_Ps(P_est,P_est1);
e01 = epipole_from_Ps(P_est1,P_est);
 pp = [0;0;1];
[Fb] = bougnoux_kruppa_f(F10,e10,e01,pp,pp);
[Fh,~] = hartley_kruppa_f(F10,e10,pp,pp);
