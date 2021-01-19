% 运行active_input之后，plot_set用于画结果图
% 运行前需移除CPLEX的路径，防止与mpt3冲突
X=Polyhedron('lb',[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],'ub',[1,1,1,1,1,1,1,1,1,1,1]);
P1=G_1*X+Y_c1;
P2=G_2*X+Y_c2;
P3=G_3*X+Y_c3;
P4=G_4*X+Y_c4;
P5=G_5*X+Y_c5;

plot(P1,'color','r')
hold on
plot(P2,'color','g')
hold on
plot(P3,'color','b')
hold on
plot(P4,'color','y')
hold on
plot(P5,'color','m')