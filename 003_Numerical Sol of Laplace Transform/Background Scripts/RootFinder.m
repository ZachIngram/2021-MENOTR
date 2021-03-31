%Root Finder
function [XbarFit,Xf]=RootFinder(Xt,Mt,Xbar)
 

 
 
% Enter binding equation here. Below, the conservation of mass equation is
% given as ImpFunc = @(Xf) Xt-Xf-(...)*Mt. In this function (...) is where
% the equation for Xbar is placed, For a single site binding
% equation this would yeild ImpFunc = @(Xf) Xt-Xf-((K*Xf)/(1+K*Xf))*Mt.
 
Xt_string = num2str(Xt);
Mt_string = num2str(Mt);
 
    
Build_Imp_Func = strcat('@(Xf) ',Xt_string,'-Xf-(',Xbar,')*',Mt_string);
 
 
ImpFunc = str2func(Build_Imp_Func);
 
 
 
 
% Solves the for the roots of the above equation to yield Xf. This is
% bounded between 0 and Xt.
Xf = fzero(ImpFunc,[0 Xt]);
 
% Enter Xbar equation here as well. Using the determined Xf, Xbar is now calculated
XbarFit=eval(Xbar);
 

