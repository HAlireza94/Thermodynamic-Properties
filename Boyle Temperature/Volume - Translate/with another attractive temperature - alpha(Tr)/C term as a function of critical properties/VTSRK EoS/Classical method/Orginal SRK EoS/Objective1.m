function [X,fval,exitflag,output]=Objective1(P,R,Tc,a,b,k,c,PARA0)

options = optimset('Display','on','MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-2,'TolFun',1e-2);
%'OutputFcn',@functn,
    
    
    [X,fval,exitflag,output]=fminsearch(@(PARA0) objfun(k,a,b,Tc,R,P,c,PARA0),PARA0,options);

        
%     function stop=functn(PARA0,optimValues)
%         stop=false;
%         if abs(optimValues.fval)<1e-05 || optimValues.fval==0
%           stop=true;
%         end
%     end



function f=objfun(k,a,b,Tc,R,P,c,PARA0)
f=(PARA0(1).^3)+(((3*((c*P)./(R*PARA0(2))))-1)*PARA0(1).^2)+(((3*(((c*P)./(R*PARA0(2))).^2))-(((b*P)./(R*PARA0(2))).^2)-(2*((c*P)./(R*PARA0(2))))-(((b*P)./(R*PARA0(2))))+((((1+(k*(1-((PARA0(2)./Tc).^0.5)))).^2)*a*P)./((R^2)*(PARA0(2).^2)))).*PARA0(1))+((((c*P)./(R*PARA0(2))).^3)-((((b*P)./(R*PARA0(2))).^2)*((c*P)./(R*PARA0(2))))-(((c*P)./(R*PARA0(2))).^2)-(((b*P)./(R*PARA0(2)))*((c*P)./(R*PARA0(2))))+(((((1+(k*(1-((PARA0(2)./Tc).^0.5)))).^2)*a*P)./((R^2)*(PARA0(2).^2)))*((c*P)./(R*PARA0(2))))-(((((1+(k*(1-((PARA0(2)./Tc).^0.5)))).^2)*a*P)./((R^2)*(PARA0(2).^2)))*((b*P)./(R*PARA0(2)))));
end

end
