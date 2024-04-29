angle = 0:0.05:2*pi
x = cos(angle)';
y = sin(angle)';
z = [x,y];
figure(1)
hold on
plot(x,y)   % circle
for i = 1:length(angle)
  plotv(z,'-')   %plot vector for each angle 
  drawnow;pause(0.2);
end