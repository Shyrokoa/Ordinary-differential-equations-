function AMB23 = AMB23(h)
close all;

Y = [0 -1];
T = -h:h:2;
global ieval
ieval = 0;

Y(1) = Y(2) - h*f(0, Y(2));
for i=2:length(T)-1
    
    % Adams-Bashforth 2rd
    Yn = Y(i) + h*f(T(i),Y(i));
    % Adams-Moulton 3th
    Y(i+1) = Y(i) + h*((5/12)*f(T(i+1),Yn)+(8/12)*f(T(i),Y(i))-(1/12)*f(T(i-1),Y(i-1)));
 
    An(i) = ya(i);
 
    % Truncation error     
    AMB(i) = trError(Y(:,i),An(:,i)); 
end
plot(T(2:end),Y(2:end), T(2:end), ya(T(2:end)));
MaxError = max(AMB)
MeanError = mean(AMB)
ieval

end

function dy = f(T,y)
global ieval;
ieval = ieval + 1;
 d=50;
dy = -d*(y - cos(T));
end

function y = ya(T)
d=50;
c1 = -1-d^2/(d^2+1);
y = c1*exp(-d*T) + d*sin(T)/(d^2+1) + d^2*cos(T)/(d^2+1);
end

function e = trError(a,b)    
e = abs(a - b); 
end 





