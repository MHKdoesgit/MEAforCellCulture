function [output, newseed] = gasdevMat( idum )

	persistent iset;
	persistent gset;

	%double fac,rsq,v1,v2;

	if( idum < 0 )
        iset = 0;
    end;
	if( iset == 0 )
        tmp = 0;
		while( tmp == 0 | rsq >= 1.0 | rsq == 0.0 )
            tmp = 1;
			[v1 newseed2] = ran1( idum );
            v1 = 2.0 * v1 - 1.0;
            idum = newseed2;
			[v2 newseed2 ] = ran1( idum );
            v2 = 2.0 * v2 - 1.0;
            idum = newseed2;
			rsq = v1*v1 + v2*v2;
        end;
            
		fac = sqrt( -2.0 * log( rsq ) / rsq );
		gset = v1 * fac;
		iset = 1;
		output = v2 * fac;
    else        
		iset = 0;
		output = gset;    
    end;
    newseed = idum;
    
