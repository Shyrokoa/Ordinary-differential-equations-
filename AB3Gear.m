function AB3Gear
Y = [0 0 -1];
err = 1e-1;
h = 0.01;

X(1) = h;
X(2) = 2*h;

global ieval
ieval = 0;

Y(2) = Y(3) - h*f(0, Y(3));
Y(1) = Y(2) - h*f(-h,Y(2));


index = 3;
step = 0;

for j = 0:h:2-(3*h)
   X(index) = X(index-1)+h;
   index = index+1;
end

for i = 3:(2/h)-1
    step = step+1;
    e = 1+err;
        while (e > err)
         % Adams-Bashforth 3rd
         Yab = Y(i)+h*((23/12)*f(X(i),Y(i))-(16/12)*f(X(i-1),Y(i-1))+(5/12)*f(X(i-2),Y(i-2)));
         % Gear 3rd
         Yg1 = (18/11)*Y(i)-(9/11)*Y(i-1)+(2/11)*Y(i-2)+(6/11)*h*f(X(i+1),Yab);
         
         % Adams-Bashforth 3rd
         Yab2 = Y(i)+0.5*h*((23/12)*f(X(i),Y(i))-(16/12)*f(X(i-1),Y(i-1))+(5/12)*f(X(i-2),Y(i-2)));
         % Gear 3rd
         Yg2 = (18/11)*Y(i)-(9/11)*Y(i-1)+(2/11)*Y(i-2)+(6/11)*0.5*h*f(X(i+1),Yab2);
         
         % Adams-Bashforth 3rd
         Yab2 = Yab2+0.5*h*((23/12)*f(X(i)+0.5*h,Yab2)-(16/12)*f(X(i-1)+0.5*h,Yab2)+(5/12)*f(X(i-2)+0.5*h,Yab2));
         % Gear 3rd
         Yg2 = Yg2+0.5*h*f(X(i+1)+0.5*h,Yab2);
                  
         err = abs(Yg1-Yg2);
         if (err > err)
             h=0.5*h;
         end
        end
        
    Y(i+1) = Yg2;
    X(i+1) = X(i)+ h;
    
    S(i) = ya(i);
    % Truncation error 
    AMB(i)=error(Y(:,i), S(:,i));
end

plot(X(3:end) ,Y(3:end));
maxError = max(AMB)
meanError = mean(AMB)

step= length(Y(3:end))
ieval

end

function dy = f(x,y)
global ieval
ieval = ieval + 1;
d = 50;
dy = -d*(y-cos(x));
end

function y = ya(x)
d = 50;
c1 = -1-d^2/(d^2+1);
y = c1*exp(-d*x)+d*sin(x)/(d^2+1)+d^2*cos(x)/(d^2+1);
end

function z = error(a,b)
    z = abs(a-b);
end