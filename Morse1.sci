 function [ECmInv,EeV] = HCl(var,J)
 //Written by Aditi and O.S.K.S. Sastri
 //var is vector consisting of De, b and Re values.
 //Typical input for var for HCl is [5,1,1.27455]
 //J is rotational quantum number. 
 //J = 0 gives pure vibrational levels
 //J = 1 gives energy eigen values corresponding to 
 //first excited rotational for v = 0, 1, 2,...
 a0 = 6;  // infinite square well width in Angstroms 
 N0 = 150; //Number of basis functions
 [V,E10] = potential(J,var); //Defines the potential function
 h = hmatrix(a0,N0,V,E10); //Determines the hmatrix 
 evals = spec(h);  //spec command used for obtaining the eigen values
 EeV = gsort(evals); //gsort command used for sorting data  in descending order
 ECmInv = EeV/(1239.84193*10^(-7));//energies in wavenumbers
 endfunction
 
 function [V,E10] = potential(J,var,a0)
 // Model parameters
 De = var(1); //Molecular dissociation energy expressed in eV
 b = var(2); //parameter in Morse potential definition expressed in Angstroms^(-1)
 Re = var(3); //equilibrium bond length in Angstroms
 R= 0.6:0.001:a0; // discretising distance parameter
 // Defining Morse potential
 v = De*((exp(-2*b*(R-Re)))-(2*exp(-b*(R-Re))));
 mu = 0.9796*931.49410*10^6; // reduced mass of molecule in eV
 hbarc = 1973.29; //value is in eV-Angstroms
 //Defining centrifugal potential for rotational term
 vcf = (J*(J+1)*hbarc^2./(2*mu*R.^2));
 V = v+vcf;
 plot(R,V); 
 // ground state energy of infinite square well potential in eV
 E10 = (%pi^2*hc^2)/(2*mu*(a0^2)); 
 endfunction
 
 function [h] = hmatrix(a0,N0,V,E10)
 h = zeros(N0,N0);
 for m = 1:N0
 h(m,m) = m^2*E10 + Vmm(m,a0,V); //Diagonal elements
    for n = m+1:N0
      h(m,n) = Vnm(n,m,a0,V); //Non-Diagonal elements
      h(n,m) = h(m,n);
    end
 end
 endfunction
 
 function [I1] = Vmm(m,a0,V)
 R = 0.6:0.001:a0;
 c1 = 1-cos(2*m*%pi*R/a0);
 f1 = c1.*V/a0;
 I1 = intsplin(R,f1); // Integration using spline interpolation
 endfunction
 
 function [I2] = Vnm(n,m,a0,V)
 R = 0.6:0.001:a0;
 c1 = cos((n-m)*%pi*R/a0);
 c2 = cos((n+m)*%pi*R/a0);
 f2 = V.*(c1-c2)/a0;
 I2 = intsplin(R,f2); //Integration using spline interpolation 
 endfunction
 
