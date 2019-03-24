function 2orderODE.m
h = 0.0001;   %[0.1,0.01,0.001,0.0001]
x(1,:) = [1]; x(2,:) = [2];
y(1,:) = [1]; y(2,:) = [2];
z(1,:) = [1]; z(2,:) = [2];
T = [0];
t = 0; i = 1;

while(t < 3)
    %Metoda Eulera jawna
    x(:,i+1) = x(:,i) + h * f(t, x(:,i));
    %Metoda Heuna
    y(:,i+1) = y(:,i) + h/2*(f(t,y(:,i))+f(t+h,y(:,i)+h*f(t, y(:,i))));       
    %Metoda Runge_Kutta_4th
    k1 = f(t,z(:,i));
    k2 = f(t+h/2,z(:,i)+h/2*k1);
    k3 = f(t+h/2,z(:,i)+h/2*k2);
    k4 = f(t+h,z(:,i)+h*k3);
    z(:,i+1) = z(:,i)+h/6*(k1+2*k2+2*k3+k4);                                   
  
    t = t+h; i = i+1;
    T(i) = t;
end
format longG
EndVal = [x(1,i), x(2,i)
          y(1,i), y(2,i)
          z(1,i), z(2,i)]

figure(1)
hold on
plot(T,x(1,:),'b-', T,x(2,:),'r-')
title('Euler explicit dla h = 0.0001;');
legend('y1','y2','Location','northwest');
hold off
figure(2)
hold on
plot(T,y(1,:),'b-', T,y(2,:),'r-')
title('Heun dla h = 0.0001;');
legend('y1','y2','Location','northwest');
hold off
figure(3)
hold on
plot(T,z(1,:),'b-', T,z(2,:),'r-')
title('Runge-Kutta 4th dla h = 0.0001;');
legend('y1','y2','Location','northwest');
hold off
end

function dy = f(t,y)
   dy = [y(2)
       -0.003*y(1)*y(1)-y(2)+5];
end
