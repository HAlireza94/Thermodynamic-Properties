function [X,fval,exitflag,output]=Objective_Nasrifar(P,R,Tc,a,b,k,Tpt,PARA0)

options = optimset('Display','on','MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-2,'TolFun',1e-2);
%'OutputFcn',@functn,
    
    
    [X,fval,exitflag,output]=fminsearch(@(PARA0) objfun(k,a,b,Tc,R,P,Tpt,PARA0),PARA0,options);

        
%     function stop=functn(PARA0,optimValues)
%         stop=false;
%         if abs(optimValues.fval)<1e-05 || optimValues.fval==0
%           stop=true;
%         end
%     end



function f=objfun(k,a,b,Tc,R,P,Tpt,PARA0)
f=((PARA0(1)).^3)-((1-((b*P)/(R.*(PARA0(2)))))*(PARA0(1)).^2)+((((a*((1+(k*(1-(((PARA0(2)-Tpt)./(Tc-Tpt)))^0.5))).^2)*P)/((R^2)*((PARA0(2)).^2)))-(2*((b*P)/(R*(PARA0(2)))))-(3*(((b^2)*(P^2))/((R^2)*((PARA0(2)).^2)))).*(PARA0(1)))+(((a*b*P*P*(1+(k*(1-((((PARA0(2)-Tpt)./(Tc-Tpt)))^0.5)))^2)/((R^3)*((PARA0(2)).^3))-((((b^2)*(P^2))/((R^2)*((PARA0(2)).^2))))-((((b^3)*(P^3))/((R^3)*((PARA0(2))^3))))))));
end

end

