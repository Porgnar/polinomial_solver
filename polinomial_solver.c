#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


/////////////////////////////////////////////////////////////////
//   Function calculates polinomial f(x) to a given x input   //
////////////////////////////////////////////////////////////////

// as an input it needs to know x
// for calculation speed reasons we make a main wide array for the various powers of x needed in calculation, we pass it to the function
// we pass the coefficients given to the program as input parameters (the number of these determines the order of our polinomial automatically)
// we pass the number of arguments of main, necessary for the program to know the order of polinomial, to know how much to iterate for each part


// our program takes inputs as coefficients of each order of x, in decreasing order, for a certain step in polyfx we will need to flip this order virtually

long double polyfx(long double x, long double xnth[], long double coeff[], int argc){
		
		long double fx=0;
		int l=argc-2;	//counter to flip order of x
		
	
		for(int d=0;d<argc-1;d++){
			
			xnth[d]=1.0;    //start each order of x as 1, this will be multiplied by our x n times, to reach x^n
		
		}
		
		for(int j=0;j<argc-2;j++){
			
			for(int k=l;k>0;k--){
				
				xnth[j] = xnth[j] * x;		//we calculate x^n, where n is an element of [ highest_order ; 0 ] integer interval
		
			}
			
			l--; //we decrease the order of x
		}
		
		
		for(int j=0;j<argc-1;j++){
			
			fx += coeff[j] * xnth[j];	//we calculate the value of f(x) by summing over each order for the x^n * coefficent product
			
		}
	
	
		return fx;
};


//////////////
//   MAIN   //
//////////////


int main(int argc, char *argv[]){
	
	// variable declarations
	
	long double coeff[argc-1];
	long double xnth[argc-1];
	
	long double lastfx = 0.0;
	long double lastx = 0.0;
	
	
	long double xinterest[2*(argc-2)];
	long double fxinterest[2*(argc-2)];
	int ic=0;
	
	
	long double roots[argc-2];
	int rc=0;
	
	// initialise array values with ones that statistically speaking should never interfere with calculations
	
	for(int i=0;i<(2*(argc-2));i++){
	
		xinterest[i] = INT_MAX;
		fxinterest[i] = INT_MAX;
		
	}
	
	// read the argv[] array for input parameters, which are our coefficients
	
	for(int i=0;i<argc-1;i++){
		
		xnth[i]=1;
		coeff[i] = atof(argv[i+1]);
		
	}
	
	
	
	
	// by adjusting the step (efficient to do it in a way of 0.5^n) we can make the initial solution interval scanning more detailed (necessary for more rapidly changing functions, or when certain roots are very close)
	for(long double i=-100.0;i<100.0;i+=0.0009765625){
			
		
		long double fx;
		
		fx = polyfx(i, &xnth[0], &coeff[0], argc);
		
		
		
		// if we find 0 explicitly by just scanning, we save the value at which we found it
		if(fx==0){
			
			// to keep the array format consistent and make our life easier later, we save the same values twice
			xinterest[ic] = i;
			fxinterest[ic] = 0.0;
			ic++;
			xinterest[ic] = i;
			fxinterest[ic] = 0.0;
			ic++;
			
		}
		
		// if we detect a change in the sign of f(x), we save the two x values between which our continuous f(x) passes zero
		if(((lastfx<0)&&(fx>0))||((lastfx>0)&&(fx<0))){
			
			xinterest[ic] = lastx;
			fxinterest[ic] = lastfx;
			ic++;
			xinterest[ic] = i;
			fxinterest[ic] = fx;
			ic++;	
			
		}
		
		//we update the previous record
		lastfx = fx;
		lastx = i;
		
		
	}
	
	
	// Actual root finding algorithm
	
	for(int i=0;i<(2*(argc-2));i++){
		// we saved the interesting values in pairs, so we only need the components of the array every second iteration
		if(((i%2)==1)){
			
			// if we didn't find roots we still have the initialisation values
			if((fxinterest[i-1]==INT_MAX) && (fxinterest[i]==INT_MAX)){
				roots[rc] = INT_MAX;
				rc++;
			}
			
			// if we found an explicit 0 we have it stored as such
			if((fxinterest[i-1]==0)&&(fxinterest[i]==0)){
				roots[rc] = xinterest[i];
				rc++;
			}
			
			// if we found a sign change which value is positive, the false position methodology changes depending on the direction
			if(fxinterest[i-1]<fxinterest[i]){
				long double tmpx=xinterest[i-1];
				long double ptmpx=INT_MAX;
				
				for(;;){
					long double cutoff = 0.0000000000005;
					long double tmpfx = polyfx(tmpx, &xnth[0], &coeff[0], argc);
					
					tmpx= ((tmpx*fxinterest[i])-(xinterest[i]*tmpfx))/(fxinterest[i]-tmpfx);
					
					if((ptmpx-tmpx)<cutoff){roots[rc]=tmpx; rc++; break;} // if we reach a cutoff point, we can stop the method, we reached the desired precision
					
					ptmpx=tmpx;
					
				}
			}
			
			// the other branch of the method
			if(fxinterest[i-1]>fxinterest[i]){
				long double tmpx=xinterest[i];
				long double ptmpx=INT_MAX;
				
				for(;;){
					long double cutoff = 0.0000000000005;
					long double tmpfx = polyfx(tmpx, &xnth[0], &coeff[0], argc);
					
					tmpx= ((xinterest[i-1]*tmpfx)-(tmpx*fxinterest[i-1]))/(tmpfx-fxinterest[i-1]);
					
					if((ptmpx-tmpx)<cutoff){roots[rc]=tmpx; rc++; break;} // if we reach a cutoff point, we can stop the method, we reached the desired precision
					
					ptmpx=tmpx;
					
				}
			}
			
		}		
	}
	
	//let's print enough lines for the user, so they know how many roots there SHOULD be, and get reminded that the ones not displayed are either not real, a real root has a higher than 1 multiplicity, or a real root was not found by our program!!!!
	for(int i=0;i<argc-2;i++){
		
		if(roots[i]==INT_MAX){
			
			printf("Not a real root / root not found / higher root multiplicity\n");
			
		}else{
			
			printf("Real root approximated to be: %Lf \n",roots[i]);
			
		}
		
	}
	
	
	
	return 0;
}