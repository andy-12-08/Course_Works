%%%% simulated annealing y = x^3 - 10*x^2 + x; [-5 10]
close all
clear all

i = [-5:0.1:10];
y = i.^3 - 10*i.^2 + i;
figure
plot(i,y,'LineWidth',2)
xlabel('X')
ylabel('y')


lb = -5;
ub = 10;

x0 = 3;

fun = @(x) x.^3-10*x.^2+x;
xmin = simulannealbnd(fun,x0,lb,ub);


k = 0;
kmax = 100000;
xold = lb + rand*(ub-lb);
step = (ub-lb)/10;
y = zeros(kmax-1,3);
while k < kmax
    yold = xold^3 - 10*xold^2 + xold;
    xnew = xold + randi([-1 1])*step;%(ub-xold)/(ub-lb)*rand(1);
    if xnew <= ub && xnew >= lb
        ynew = xnew^3 - 10*xnew^2 + xnew;
        diff = ynew - yold;
        
        if diff < 0
            xold = xnew;
            k = k+1;
        elseif exp(-diff/kmax)>rand
            xold = xnew;
        else 
            k = k+1;
        end
    else 
        k = k+ 1;

    y(k,1)=yold;
 y(k,2)=ynew;
 y(k,3)=diff;
    end

 
 
end
xmin
xnew
ynew
 


