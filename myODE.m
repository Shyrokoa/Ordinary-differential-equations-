function [] = myODE(d,h)
yE(1) = d; 
yH(1) = d;
yRK(1) = d;
T = 0;
t = 0; 
i = 1;
a = 1/d-1;
yL(1) = yH(1);

EulerError = 0; 
HeunError = 0; 
RKError = 0;

while(t < 2/d)
    % Lambert
    yL(i+1) = 1/(lambertw(a*exp(a-t))+1);
    % Euler explicit method
    yE(:,i+1) = yE(:,i) + h * f(t, yE(:,i));
    % Heun's method
    yH(:,i+1) = yH(:,i) + h/2*(f(t,yH(:,i))+f(t+h,yH(:,i)+h*f(t, yH(:,i))));       
    % Runge-Kutta 4th range method
    k1 = f(t,yRK(:,i));
    k2 = f(t+h/2,yRK(:,i)+h/2*k1);
    k3 = f(t+h/2,yRK(:,i)+h/2*k2);
    k4 = f(t+h,yRK(:,i)+h*k3);
    yRK(:,i+1) = yRK(:,i)+h/6*(k1+2*k2+2*k3+k4);                                   
    % Truncation error
    EulerError(i) = trError(yL(:,i),yE(:,i));
    HeunError(i) = trError(yL(:,i),yH(:,i));
    RKError(i) = trError(yL(:,i),yRK(:,i));
    
    t = t+h; 
    i = i+1;
    T(i) = t;
end


MaxError = [max(EulerError), max(HeunError), max(RKError)]
MeanError= [mean(EulerError), mean(HeunError), mean(RKError)]

plot(T,yL(1,:),'r', T,yE(1,:),'g', T,yH(1,:),'b', T,yRK(1,:),'m');
hold on
legend('Lambert','Euler explicit','Heun`s method','Runge-Kutta 4th order','Location','best');
end

function dy = f(~,y)
   dy = y^2-y^3;
end

function e = trError(yL,y)
   e = abs(yL - y);
end
