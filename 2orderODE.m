function 2orderODE.m
tFinal = 3;
yE(1,:) = 1; yE(2,:) = 2;
yH(1,:) = 1; yH(2,:) = 2;
yRK(1,:) = 1; yRK(2,:) = 2;
T = 0;
t = 0; i = 1;

while(t < tFinal)
    % Euler explicit method
    yE(:,i+1) = yE(:,i) + h * f(t, yE(:,i));
    % Heun's method
    yH(:,i+1) = yH(:,i) + h/2*(f(t,yH(:,i))+f(t+h,yH(:,i)+h*f(t, yH(:,i))));       
    % Runge-Kutta 4th range method
    K1 = f(t,yRK(:,i));
    K2 = f(t+h/2,yRK(:,i)+h/2*K1);
    K3 = f(t+h/2,yRK(:,i)+h/2*K2);
    K4 = f(t+h,yRK(:,i)+h*K3);
    yRK(:,i+1) = yRK(:,i)+h/6*(K1+2*K2+2*K3+K4);                                   
  
    t = t+h; 
    i = i+1;
    T(i) = t;
end


% y1(3) and y2(3)
finishValue = [ yE(1,i), yE(2,i)
                yH(1,i), yH(2,i)
                yRK(1,i), yRK(2,i)];

figure(1)
hold on
plot(T,yE(1,:),'g-', T,yE(2,:),'r-')
title('Euler explicit method for h = 0.0001;');
legend('y1','y2','Location','best');
grid on
hold off

figure(2)
hold on
plot(T,yH(1,:),'g-', T,yH(2,:),'r-')
title('Heun`s method for h = 0.0001;');
legend('y1','y2','Location','best');
grid on
hold off

figure(3)
hold on
plot(T,yRK(1,:),'g-', T,yRK(2,:),'r-')
title('Runge-Kutta 4th range method for h = 0.0001;');
legend('y1','y2','Location','best');
grid on
hold off
end

function dy = f(~,y)
   dy = [y(2) ; -0.003*(y(1)^2)-y(2)+5];
end
