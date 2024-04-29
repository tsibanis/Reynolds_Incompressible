t1 = readmatrix('t_30_30.dat'); 
t5 =readmatrix('t_90_30.dat'); 
t2 = readmatrix('t_100_100.dat'); 
t3 = readmatrix('t_300_100.dat'); 
t4 = readmatrix('t_300_300.dat'); 
t6=readmatrix('t_500_500.dat'); 
%xlim([0 7])

ylim([-100 1.4E7])
plot(t1(:,1),t1(:,2),'LineWidth',1)
hold on
%plot(t5(:,1),t5(:,2))
%hold on
plot(t2(:,1),t2(:,2),'LineWidth',1)
hold on
%plot(t3(:,1),t3(:,2))
%hold on
plot(t4(:,1),t4(:,2),'LineWidth',1)
hold on
plot(t6(:,1),t6(:,2),'LineWidth',1)

legend('30x30','100x100','300x300','500x500')'
%legend('Circumferential Pressure Distribution for z=L/2 ')
xlabel('Theta (rad)')
ylabel('Pressure (Pa)')