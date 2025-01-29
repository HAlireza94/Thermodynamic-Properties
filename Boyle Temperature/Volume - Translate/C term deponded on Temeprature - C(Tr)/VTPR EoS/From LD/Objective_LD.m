function [X,fval,exitflag,output]=Objective_LD(P,R,Tc,a,b,k,Cc,Beta,Gamma,PARA0)

options = optimset('Display','on','MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-2,'TolFun',1e-2);
%'OutputFcn',@functn,
    
    
    [X,fval,exitflag,output]=fminsearch(@(PARA0) objfun(k,a,b,Tc,R,P,Cc,Beta,Gamma,PARA0),PARA0,options);

        
%     function stop=functn(PARA0,optimValues)
%         stop=false;
%         if abs(optimValues.fval)<1e-05 || optimValues.fval==0
%           stop=true;
%         end
%     end



function f=objfun(k,a,b,Tc,R,P,Cc,Beta,Gamma,PARA0)
f=(PARA0(1)^3)+((((b*P)/(R*PARA0(2)))-1+(3*(((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2)))))*PARA0(1)^2)+(((((a*P*((1+(k*(1-((PARA0(2)/Tc)^0.5))))^2))/((R^2)*(PARA0(2)^2)))-(2*((b*P)/(R*PARA0(2))))-(3*(((b*P)/(R*PARA0(2)))^2))+(3*((((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2)))^2))+(2*((b*P)/(R*PARA0(2)))*(((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2))))-(2*(((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2)))))*PARA0(1)))+(((((b*P)/(R*PARA0(2)))^3)+(((b*P)/(R*PARA0(2)))^2)-(((a*P*((1+(k*(1-((PARA0(2)/Tc)^0.5))))^2))/((R^2)*(PARA0(2)^2)))*((b*P)/(R*PARA0(2))))+(((b*P)/(R*PARA0(2)))*((((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2)))^2))+((((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2)))^3)-(3*(((b*P)/(R*PARA0(2)))^2)*(((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2))))-((((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2)))^2)-(2*((b*P)/(R*PARA0(2)))*(((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2))))+(((a*P*((1+(k*(1-((PARA0(2)/Tc)^0.5))))^2))/((R^2)*(PARA0(2)^2)))*(((Cc*(Beta+((1-Beta)*exp(Gamma*(abs(1-(PARA0(2)/Tc)))))))*P)/(R*PARA0(2))))));
end

end
