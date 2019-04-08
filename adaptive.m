function adaptive = adaptive(stiff)
global ieval;
global d;
d = stiff;
ieval = 0;
Y=[0 -1];
err=1e-3;
h=0.001;

Y(1)=Y(2)-h*f(0,Y(2));
X(1)=h;

index=2;

for j=0:h:2-(2*h)
   X(index)=X(index-1)+h;
   index=index+1;
end

for i=2:(2/h)-1
    e = 1+err;
        while (e > err)
            
         %  Adams-Bashforth 2rd
         Yab=Y(i)+h*(3/2*f(X(i),Y(i))-1/2*f(X(i-1),Y(i-1)));
         %  Adams-Moulton 3th
         Yam=Y(i)+h*(5/12*f(X(i+1),Yab)+2/3*f(X(i),Y(i))-1/12*f(X(i-1),Y(i-1)));
         
         %  Adams-Bashforth 2rd
         Yab2=Y(i)+0.5*h*(3/2*f(X(i),Y(i))-1/2*f(X(i-1),Y(i-1)));
         %  Adams-Moulton 3th
         Yam2=Y(i)+0.5*h*(5/12*f(X(i+1),Yab2)+2/3*f(X(i),Y(i))-1/12*f(X(i-1),Y(i-1)));
         
         %  Adams-Bashforth 2rd
         Yab2=Yab2+0.5*h*(3/2*f(X(i)+0.5*h,Yab2)-1/2*f(X(i-1)+0.5*h,Yab2));
         %  Adams-Moulton 3th
         Yam2=Yam2+0.5*h*(5/12*f(X(i+1)+0.5*h,Yab2)+2/3*f(X(i)+0.5*h,Yam2)-1/12*f(X(i-1)+0.5*h,Yam2));
         
         e=abs(Yam2-Yam);
         if (e > err)
             h=0.5*h;
         end
        end
      
        Y(i+1)= Yam2;
        X(i+1)=X(i)+ h;
        
        S(i) = ya(i);
        % Truncation error     
        AMB(i) = trError(Y(:,i),S(:,i)); 
end

MaxError = [max(AMB)] 
MeanError= [mean(AMB)]
ieval

plot(X(2:end) ,Y(2:end) ,X(2:end), ya(X(2:end)));
end


function e = trError(a,b)    
e = abs(a - b); 
end 

function dy= f(x,y)
global ieval;
global d;
ieval = ieval + 1;
dy=-d*(y-cos(x));
end

function y=ya(x)
global d;
c1=-1-d^2/(d^2+1);
y=c1*exp(-d*x)+d*sin(x)/(d^2+1)+d^2*cos(x)/(d^2+1);
end
