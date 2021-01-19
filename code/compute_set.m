function [G_psi,Y_center,N] = compute_set(A,B)

C1=[1 0];
Bw1=[1 0; 0 0];
Dv1=1;
% u0 = sym('u0',[2 1]);
% u1 = sym('u1',[2 1]);
u0=[0;0];
u1=[0;0];
Gx=[0.2 0;0 0.2];
Gw=[0.5 0;0 0.5];
Gv=0.2;
Cx=[-3;-3];
Cw=[0;0];

A_tilde=[eye(2);A;A*A];
B_tilde=[zeros(2,4);B zeros(2,2);A*B B];
Bw_tilde=blkdiag(zeros(2),Bw1,A*Bw1+Bw1);
C_tilde=blkdiag(C1,C1,C1);
Gw_tilde=blkdiag(Gw,Gw,Gw);
G_phi=[A_tilde*Gx Bw_tilde*Gw_tilde];
Dv_tilde=blkdiag(Dv1,Dv1,Dv1);
Gv_tilde=blkdiag(Gv,Gv,Gv);
G_psi=[C_tilde*G_phi Dv_tilde*Gv_tilde];

mid=sum((abs(G_psi)),2);
G_psi=diag(mid);

X_c0=Cx;
X_c1=A*X_c0+B*u0;
X_c2=A*X_c1+B*u1;
X_center=[X_c0;X_c1;X_c2];

Y_c0=C1*X_c0;
Y_c1=C1*X_c1;
Y_c2=C1*X_c2;
Y_center=[Y_c0;Y_c1;Y_c2];

N=C_tilde*B_tilde;
end

