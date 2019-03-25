syms a b T0 x y p H C  g x0 F0 T0 h0 v p2  rgg gg gr lgg lgt M mgg mgt rgt  th1 th2 th3 th4 th5 tha th6 T01 T12 T23 T34 T45 T56 f Fl y0 z(x) x1 x2
M=1000;g=9.8;mgg=10;mgt=100;p2=1025;rgg=0.025;lgg=1;rgt=0.15;lgt=1;p=7;
mq=1200;
%v=24;
%h0=0.74883;
 %v=24,h0=0.74883;    %反复设置h0,使答案逼近深度18m % h0=0.748828125;
% v=12; h0= 0.7347;
v=36 ; h0=0.7698
 
F0=1.25*v^2 *(2-h0); %风力
Ffl = p2*g*pi*h0; %3.1556*10^4*h   %浮力
th1 = atan((1.25*v^2 *(2-h0))/(3.1556*10^4*h0 - M*g)); %thi 第i节钢管的夹角
th2 = atan(F0/(Ffl-M*g + (p2*g*pi*rgg^2*lgg - mgg*g)));
th3 = atan(F0/(Ffl-M*g + 2*(p2*g*pi*rgg^2*lgg - mgg*g)));
th4 = atan(F0/(Ffl-M*g + 3*(p2*g*pi*rgg^2*lgg - mgg*g)));
T45 = F0 / sin(th4);    
th5 = atan(F0/(Ffl-M*g + 4*(p2*g*pi*rgg^2*lgg - mgg*g) + p2*g*pi*rgt^2*lgt - mgt*g - mq*g)); %锚链上端切线方向
th7= atan(F0/(Ffl-M*g + 4*(p2*g*pi*rgg^2*lgg - mgg*g) + p2*g*pi*rgt^2*lgt - mgt*g - mq*g - 22.05*7)); %无用
T56 = F0 / sin(th5);

ap = F0/(p*g);
th6 = solve('T56*sin(th5-th6) = T45*sin(th6-th4) + mq*g*sin(th6)',th6); %钢桶实际倾角 力矩
th6=th6(1);
%th6 =eval(th6(1));


f(x)=ap*(cosh(x/ap) - 1); %悬链线方程

%z(x1)=vpa(subs(diff(f(x)),x,x1));

% x1 = solve('atan(sinh(0.28614093483828451376137126812156*x1))==pi/2-th5',x1) %v=12时
 %x1=solve('sinh(x1/ap)==tan(pi/2-th5)',x1);
% x1=asinh(tan(pi/2 - th5))*ap
 %x1=eval(x1(1))
 x1=asinh(tan(pi/2 - th5))*ap
 
%x2=double(solve(ap*(sinh(x1/ap)-sinh(x2/ap))-22.05,x2));
x2=solve('ap*(sinh(x1/ap)-sinh(x2/ap))==22.05',x2);%根据锚链长度求锚链下端坐标x2
x2=eval(x2(1))

%钢桶方程
Y1(x)=tan(pi/2-th6)*(x-x1)+f(x1)
%第四节钢管方程
Y2(x)=tan(pi/2-th4)*(x-x1-sin(th6))+Y1(x1+sin(th6))
%第三节钢管方程
Y3(x)=tan(pi/2-th3)*(x-x1-sin(th6)-sin(th4))+f(x1)+cos(th6)+cos(th4)
%第二节钢管方程
Y2(x)=tan(pi/2-th2)*(x-x1-sin(th6)-sin(th4)-sin(th3))+f(x1)+cos(th6)+cos(th4)+cos(th3)
%第一节钢管方程
Y1(x)=tan(pi/2-th1)*(x-x1-sin(th6)-sin(th4)-sin(th3)-sin(th2))+f(x1)+cos(th6)+cos(th4)+cos(th3)+cos(th2)


if x2<0  %锚链松弛
    l = eval(f(x1))
      b=0
      L=22.05-ap*(sinh(x1/ap)) %着地长度
      
else
    l = eval(f(x1)-f(x2))
    b=sinh(x2/ap)
end

f2 =eval( l + lgg*(cos(th1)+cos(th2)+cos(th3)+cos(th4))+lgt*cos(th6) + h0 - 18) %计算深度误差

th6 =eval(th6(1));
h= l + lgg*(cos(th1)+cos(th2)+cos(th3)+cos(th4))+lgt*cos(th6) + h0 ;        %预测海深
r=x1 + lgg*(sin(th1)+sin(th2)+sin(th3)+sin(th4))+lgt*sin(th6) %浮动半径

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

axis([-10 30 0 20]) 
hold off




