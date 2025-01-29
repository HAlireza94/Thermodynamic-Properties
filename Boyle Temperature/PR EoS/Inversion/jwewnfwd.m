clc
clear
close all

    Tc=126.2;Pc=34;w=0.039;Zc=0.294;Vc=89.21;
  R=83.14472;

         k=0.37464+(1.54226*w)-(0.26992*w^2);

          a=(0.45724*(R^2)*(Tc^2))/Pc;
  b=(0.07780*R*Tc)/Pc;

  oma=0.45724;omb=0.0778;
         
  
  T=linspace(0.7*Tc,5.05*Tc,100)';
  
  p=Pc;
  
  X=zeros(100,1);
  Tr=zeros(100,1);
  Pr=zeros(100,1);
  for i=1:100
      
      Tr(i)=T(i)/Tc;
      
      v(i)=((R*T(i))/p);
      x0(i)=(v(i)/Vc);
      
        alpha(i)=((1+(k*(1-(Tr(i)^0.5))))^2);
        first_order_alpha(i)=(-k/sqrt(Tr(i)))*sqrt(alpha(i));
        
        F=@(x) (Tr(i)*((1/((x*Zc)-omb))-((oma*first_order_alpha(i))/(((x*Zc)^2)+(2*x*Zc*omb)-omb^2))))+...
            (x*(((-Tr(i)*Zc)/(((x*Zc)-omb)^2))+((2*oma*alpha(i)*Zc*(x+omb))/((((x*Zc)^2)+(2*omb*Zc*x)-b^2)^2))));
        
        X(i)=fzero(F,x0(i));
        
        
        Pr(i)=(Tr(i)/((X(i)*Zc)-omb))-((oma*alpha(i))/(((X(i)*Zc)^2)+(2*omb*X(i)*Zc)-b^2));
        
                                                                        
        
  end
        