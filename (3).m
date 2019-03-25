syms a b T0 p H C x y g x0 F0 T0 h0 v p2  rgg gg gr lgg lgt M mgg mgt rgt  th1 th2 th3  tha  T01 T12 T23 T34  f Fl y0 z(x) x1  f2
M=1000;g=9.8;mgg=10;mgt=100;p2=1025;rgg=0.025;lgg=1;rgt=0.15;lgt=1;p=19.5;mq=2000;v=36;v2=1.5

istrue=0;

for p=[ 3.2 7 12.5 19.5 28.12]
    for h0=0.6:0.01:2

     syms th4 th5 th6 T45 T56 x2
    
F0=1.25*v^2 *(2-h0)+374*v2^2*h0*2; %风力

Ffl = p2*g*pi*h0; %3.1556*10^4*h   %浮力
th1 = atan((1.25*v^2 *(2-h0))/(3.1556*10^4*h0 - M*g));  %thi 第i节钢管的夹角
th2 = atan(F0/(Ffl-M*g + (p2*g*pi*rgg^2*lgg - mgg*g)));
th3 = atan(F0/(Ffl-M*g + 2*(p2*g*pi*rgg^2*lgg - mgg*g)));
th4 = atan(F0/(Ffl-M*g + 3*(p2*g*pi*rgg^2*lgg - mgg*g)));
T45 = F0 / sin(th4);    
th5 = atan(F0/(Ffl-M*g + 4*(p2*g*pi*rgg^2*lgg - mgg*g) + p2*g*pi*rgt^2*lgt - mgt*g - mq*g)); %钢桶倾角
th7= atan(F0/(Ffl-M*g + 4*(p2*g*pi*rgg^2*lgg - mgg*g) + p2*g*pi*rgt^2*lgt - mgt*g - mq*g - 22.05*7)); %无用
T56 = F0 / sin(th5);

ap = F0/(p*g);

th6 = solve('T56*sin(th5-th6) = T45*sin(th6-th4) + mq*g*sin(th6)',th6); %钢桶实际倾角 力矩
th6=th6(1);
%th6 =eval(th6(1));


f(x)=ap*(cosh(x/ap) - 1); %悬链线方程

x1=asinh(tan(pi/2 - th5))*ap
%z(x1)=vpa(subs(diff(f(x)),x,x1));
 %x1 = solve('atan(sinh(x1/ap))==pi/2 - th5',x1); %求锚链上端坐标x1
 %x1=eval(x1(1))
 
%x2=double(solve(ap*(tan(pi/2 - th5)-sinh(x2/ap))-22.05,x2));
x2=solve('ap*(sinh(x1/ap)-sinh(x2/ap))==32',x2);%根据锚链长度求锚链下端坐标x2
x2=eval(x2(1))


if x2<0 %锚链松弛
    l = eval(f(x1))
    b=0
    L=32-ap*(sinh(x1/ap)) %着地长度
else
    l = eval(f(x1)-f(x2))
    b=sinh(x2/ap) %锚点与海床的夹角
end


f2 =abs(eval(l + lgg*(cos(th1)+cos(th2)+cos(th3)+cos(th4))+lgt*cos(th6(1)) + h0 - 18)/18) %计算深度误差

th6 =eval(th6(1))

if f2<0.05 &&th5>0 
    
   h= l + lgg*(cos(th1)+cos(th2)+cos(th3)+cos(th4))+lgt*cos(th6) + h0 ;        %预测海深
   r=abs(x1) + lgg*(sin(th1)+sin(th2)+sin(th3)+sin(th4))+lgt*sin(th6) %浮动半径
   istrue=1
   
   break;
end

    end
    
   if x2<0
    
 fplot(@(x)0,[-L 0],'b'); 
 hold on
 fplot(@(x)ap*(cosh(x/ap) - 1),[0 x1],'r'); 
 
else
  fplot(@(x)ap*(cosh(x/ap) - 1),[0 x2],'g'); 
  hold on
   fplot(@(x)ap*(cosh(x/ap) - 1),[x2 x1],'r'); 
    
end

 hold on
 fplot(@(x)(tan(pi/2-th6)*(x-x1)+f(x1)),[x1 x1+sin(th6)],'g'); 
 hold on
 fplot(@(x)tan(pi/2-th4)*(x-x1-sin(th6))+f(x1)+cos(th6),[x1+sin(th6) x1+sin(th6)+sin(th4)],'r'); 
 hold on
 fplot(@(x)tan(pi/2-th3)*(x-x1-sin(th6)-sin(th4))+f(x1)+cos(th6)+cos(th4),[x1+sin(th6)+sin(th4) x1+sin(th6)+sin(th4)+sin(th3)],'g'); 
 hold on
 fplot(@(x)tan(pi/2-th2)*(x-x1-sin(th6)-sin(th4)-sin(th3))+f(x1)+cos(th6)+cos(th4)+cos(th3),[x1+sin(th6)+sin(th4)+sin(th3) x1+sin(th6)+sin(th4)+sin(th3)+sin(th2)],'r'); 
 hold on
 fplot(@(x)tan(pi/2-th1)*(x-x1-sin(th6)-sin(th4)-sin(th3)-sin(th2))+f(x1)+cos(th6)+cos(th4)+cos(th3)+cos(th2),[x1+sin(th6)+sin(th4)+sin(th3)+sin(th2) x1+sin(th6)+sin(th4)+sin(th3)+sin(th2)+sin(th1)],'g'); 
 hold on

end
 axis([-10 30 0 20]) 
     hold off
   


