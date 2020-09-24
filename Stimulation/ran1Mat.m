function [output, newseed] = ran1Mat( idum )

    persistent iy;
    persistent iv;
    
    IA  = 16807;
    IM = 2147483647;
    AM = (1.0/IM);
    IQ = 127773;
    IR = 2836;
    NTAB = 32;
    NDIV = (1+floor((IM-1)/NTAB));
    EPS = 1.2e-7;
    RNMX = (1.0-EPS);

	%int j;
	%long k;
	%static long iy=0;
	%static long iv[NTAB];
	%double temp;

    if( idum <= 0 )
        iy = 0;
        iv = zeros( 1, NTAB );
    end;
	if( idum <= 0 | iy == 0 )
		if( -1*idum < 1 )
            idum = 1;
        else
            idum = -1 * idum;
        end;
		for j = NTAB + 7 : -1: 0;
			k = floor( idum / IQ );
			idum = IA * ( idum - k*IQ ) - IR*k;
			if( idum < 0)
                idum = idum + IM;
            end;
			if( j < NTAB )
                iv(j+1) = idum;
            end;
        end;
		iy = iv(1);
    end;
    
	k = floor( idum/IQ );
	idum = IA * ( idum - k*IQ ) - IR*k;
	if( idum < 0)
        idum = idum + IM;
    end;
	j = floor( iy/NDIV );
	iy = iv(j+1);
	iv(j+1) = idum;
    temp = AM * iy;
    newseed = idum;
	if( temp > RNMX )
        output = RNMX;
    else
        output = temp;
    end;

    