syms a b T0 p H C x y g x0 F0 T0 h0 v p2  rgg gg gr lgg lgt M mgg mgt rgt  th1 th2 th3  tha  T01 T12 T23 T34  f Fl y0 z(x) x1  f2  s mq mq0
M=1000;g=9.8;mgg=10;mgt=100;p2=1025;rgg=0.025;lgg=1;rgt=0.15;lgt=1; p=28.12
v=36;v2=1.5
%s=22;mq=2000
s=26;mq=2400

istrue=0;

     for h0=0.8:0.01:1.1
       % for s=20:0.1:30
            
    
     syms th4 th5 th6 T45 T56 x2 
    
F0=1.25*v^2 *(2-h0); %风力

Ffl = p2*g*pi*h0; %3.1556*10^4*h   %浮力
th1 = atan((1.25*v^2 *(2-h0))/(3.1556*10^4*h0 - M*g));  %thi 第i节钢管的夹角
th2 = atan(F0/(Ffl-M*g + (p2*g*pi*rgg^2*lgg - mgg*g)));
th3 = atan(F0/(Ffl-M*g + 2*(p2*g*pi*rgg^2*lgg - mgg*g)));
th4 = atan(F0/(Ffl-M*g + 3*(p2*g*pi*rgg^2*lgg - mgg*g)));
T45 = F0 / sin(th4); 
ap = F0/(p*g);



th5 = atan((F0+374*v2^2*cos(th4)*0.3*1)/(Ffl-M*g + 4*(p2*g*pi*rgg^2*lgg - mgg*g) + p2*g*pi*rgt^2*lgt - mgt*g - mq*g)); %
T56 = F0 / sin(th5);
th6 = solve('T56*sin(th5-th6) = T45*sin(th6-th4) + mq*g*sin(th6)',th6); %钢桶实际倾角 力矩
th6=th6(1);
%th6 =eval(th6(1));

f(x)=ap*(cosh(x/ap) - 1); %悬链线方程

x1=asinh(tan(pi/2 - th5))*ap

x2=solve('ap*(sinh(x1/ap)-sinh(x2/ap))==s',x2);%根据锚链长度求锚链下端坐标x2
x2=eval(x2(1))


if x2<0 %锚链松弛
    l = eval(f(x1))
    b=0
    L=s-ap*(sinh(x1/ap)) %着地长度
else
    l = eval(f(x1)-f(x2))
    b=sinh(x2/ap) %锚点与海床的夹角
    L=0;
end

f2 =abs(eval((l + lgg*(cos(th1)+cos(th2)+cos(th3)+cos(th4))+lgt*cos(th6(1)) + h0 - 18))/18) %计算深度误差

th6 =eval(th6(1))

 h= l + lgg*(cos(th1)+cos(th2)+cos(th3)+cos(th4))+lgt*cos(th6) + h0 ;        %预测海深
   r=abs(x1) + lgg*(sin(th1)+sin(th2)+sin(th3)+sin(th4))+lgt*sin(th6) %浮动半径
   
if f2<0.05 &&th5>0 && rad2deg(b)<=16 && rad2deg(th6)<=5
   istrue=1
   break;
   
end
     

            end
%      if istrue==1
%          break;
%      end
%         end
        
%        if istrue==1
%           break;
%       end
%         
%         
%      end
%      
  
   
  
   
   
%     if istrue==1
%     break;
%    end
 % end
 
  %T=1000;minT=1;num=1000;xishu=0.95;
  
%H(mq)=atan(2*F0/((18081*pi)/80 - (49*mq)/5 + 10045*pi*h0 - 195151997835498253/17592186044416+10045*pi*h0 - 88267298721398565/8796093022208+mq*g))

    %模拟退火算法
%    %mq0=1200
%     %while T>minT
%         
%         dE=1000*(H(mq0+20)-H(mq0)) 因为精度问题，所以差别较小时dE值被取为0，所以扩大dE
%         
%         if dE<=0
%             mq0=mq0+20;
%         else
%             if exp(-dE/T)>random('Normal',0,1)
%                 mq0=mq0+20
%             end
%         end
%         
%         T=T*xishu;
%         
%         end

    
  



