clear all;
clc

A1=[0.6 0.2;-0.4 -0.2];
B1=[1 0;0 1]*1.1;
C1=[1 0];
Bw1=[1 0; 0 0];
Dv1=[1];

B2=[1 0;0 0];
B3=[0 0;0 1];
A4=[1.2 0.2;-0.4 -0.2];
A5=[2 0.2;-0.4 -0.7];

%%
u0 = sym('u0',[2 1]);
u1 = sym('u1',[2 1]);
Gx=[0.2 0;0 0.2];
Gw=[0.5 0;0 0.5];
Gv=0.2;
Cx=[-3;-3];
Cw=[0;0];
%% model 1 2 3 4 5
[G1_psi,Y1_center,N1]=compute_set(A1,B1);
[G2_psi,Y2_center,N2]=compute_set(A1,B2);
[G3_psi,Y3_center,N3]=compute_set(A1,B3);
[G4_psi,Y4_center,N4]=compute_set(A4,B1);
[G5_psi,Y5_center,N5]=compute_set(A5,B1);

%%
N=cell(1,10);
N{1}=N2-N1;
N{2}=N3-N1;
N{3}=N4-N1;
N{4}=N5-N1;
N{5}=N3-N2;
N{6}=N4-N2;
N{7}=N5-N2;
N{8}=N4-N3;
N{9}=N5-N3;
N{10}=N5-N4;
%%
G=cell(1,10);
G{1}=[G1_psi G2_psi];
G{2}=[G1_psi G3_psi];
G{3}=[G1_psi G4_psi];
G{4}=[G1_psi G5_psi];
G{5}=[G2_psi G3_psi];
G{6}=[G2_psi G4_psi];
G{7}=[G2_psi G5_psi];
G{8}=[G3_psi G4_psi];
G{9}=[G3_psi G5_psi];
G{10}=[G4_psi G5_psi];
%%
Center=cell(1,10);
Center{1}=Y1_center-Y2_center;
Center{2}=Y1_center-Y3_center;
Center{3}=Y1_center-Y4_center;
Center{4}=Y1_center-Y5_center;
Center{5}=Y2_center-Y3_center;
Center{6}=Y2_center-Y4_center;
Center{7}=Y2_center-Y5_center;
Center{8}=Y3_center-Y4_center;
Center{9}=Y3_center-Y5_center;
Center{10}=Y4_center-Y5_center;
%%
det_m=30;
tic;
p1 = binvar(10,6);
p2 = binvar(10,6);
u = sdpvar(4,1);
ksi=sdpvar(6,10);
det=sdpvar(10);
lambda=sdpvar(3,10);
mu1=sdpvar(6,10);
mu2=sdpvar(6,10);
C=[norm(u,'inf')<=9];
for q=1:10
        %添加约束
    C=[ C
        N{q}*u==G{q}*ksi(:,q)+Center{q}
        norm(ksi(:,q),'inf')<=1+det(q)
        G{q}'*lambda(:,q)==mu1(:,q)-mu2(:,q)
        1==sum(mu1(:,q)+mu2(:,q))
        0.01<=det(q)<=det_m
    ];
    for i=1:6
        C=[C;mu1(i,q)>=0;mu2(i,q)>=0];  
        C=[C;mu1(i,q)<=p1(q,i);mu2(i,q)<=p2(q,i)];
        C=[C;(-2*(1+det_m)*(1-p1(q,i)))<=(ksi(i,q)-1-det(q))<=0];
        C=[C;0<=(ksi(i,q)+1+det(q))<=(2*(1+det_m)*(1-p2(q,i)))]; 
    end
end

 %目标函数
z=norm(u,2)^2;
ops=sdpsettings('solver', 'cplex');%选择求解器
result=solvesdp(C,z,ops);%求解
 
if result.problem == 0
    uresult = value(u)
    zresult = value(z)
else
    disp('无解')
end
toc;

[G_1,Y_c1] = final_set(A1,B1,uresult);
[G_2,Y_c2] = final_set(A1,B2,uresult);
[G_3,Y_c3] = final_set(A1,B3,uresult);
[G_4,Y_c4] = final_set(A4,B1,uresult);
[G_5,Y_c5] = final_set(A5,B1,uresult);







