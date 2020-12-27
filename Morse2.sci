//Scilab code to determine model parameters for HCl molecule using VMC}
function [var_n] = VMCHCl(var_o,NumOfItns,I)
 //Written by Aditi and O.S.K.S. Sastri
 // var_o is an array consisting of model parameters given as input 
 // NumOfItns denotes number of iterations
 // I defines limits of interval range [-I,I] in which random value r is generated
 // for varying any parameter as var(i) = var(i) + r
 
 FreqExpt = [2886;5668;8347]; // Exp. vibrational frequencies in wavenumbers
 [Eno] = HCl(var_o,J); // Determining energy eigen values using matrix method
 N = length(Eno); // N stores length of eigen values vector
 // Determining first three vibrational frequencies 
 for i = 1:3
   FreqSimOld(i) = (Eno(N-i) - Eno(N));
 end
 SqrErrorOld = (FreqExpt - FreqSimOld).^2;
 chisqrold = mean(SqrErrorOld); // calculate old chi square value
 disp(chisqrold); // displaying old chi-square value
 n = 1; //initialising variable for iterative minimization
 while i <= NumOfItns 
 k = 1; // Choosing kth model parameter 
 while k <= 3
 var_n = var_o; // making old parameters as new ones
 // Monte-Carlo step
 r = (-I+2*I*rand());//random value in interval [-I,I]
 var_n(k) = var_o(k) + r; 
 [Enn] = HCl(var_n,J); // Redetermining energy eigen values using matrix method
 N = length(Enn);
 for j = 1:3
  FreqSimNew(j) = (Enn(N-j)-Enn(N)); //Redetermining vibrational freqs with new constants
 end
 SqrErrorNew = (FreqExpt - FreqSimNew).^2;
 chisqrnew = mean(SqrErrorNew); // Calculating new chi-square value 
 //Variational step
 if (chisqrnew) < (chisqrold) then 
     var_o = var_n; // Updating old variables with new ones
     chisqrmin = chisqrnew; // making new chi-square value as minimum value 
     chisqrold = chisqrnew; // updating chisqrold as chosqrnew
 else
     chisqrmin = chisqrold; //retaining old chisqr value as minimum value 
 end
 k = k+1;
 disp([FreqSimNew,FreqExpt],chisqrmin); 
 end
 n = n+1;
 end
 endfunction
