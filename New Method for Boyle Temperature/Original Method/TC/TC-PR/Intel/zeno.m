function Tr=zeno(x,RhoR,b,c,a,R,Tc,Vc)



Tr=(((R*x*Tc*RhoR)/(Vc-(b*RhoR)))...
                        -((a*((x^(N*(M-1)))*(exp(L*(1-(x^(M*N)))))))/...
                        ((((Vc/RhoR)+c)*((Vc/RhoR)+b+(2*c)))+((b+c)*((Vc/RhoR)-b))))...
                        -((Z*R*x*Tc*RhoR)/Vc));



end