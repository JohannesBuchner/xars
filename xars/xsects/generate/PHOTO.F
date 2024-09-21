c  Routine to return the photoelectric absorption cross section in cm**2
c  at energy keV for element Z.

      FUNCTION photo(keV1, keV2, Z, versn, status)

      REAL photo, keV1, keV2
      INTEGER Z, versn, status

c  Cross-section data from Henke etal, Atomic Data and Nuclear Data Tables
c  vol 27, no 1, 1982. Fits mainly by Monika Balucinska-Church and Dan McCammon
c  "Photoelectric Absorption Cross Sections with Variable Abunances"
c  Ap.J. 400, 699 (1992)

c  kaa   4/16/92
c       11/26/99     added version number so that new and old versions of
c                    the He cross-section can be accomodated.

c  Arguments :
c	keV1	r	i: Lower energy of bin in keV.
c	keV2	r	i: Upper energy of bin in keV.
c       Z       i	i: Atomic number of element
c       versn   i       i: 2 == old Marr & West He x-section
c                          3 == new Yan et al. x-section
c	status	i	r: 0 = OK
c			   1 = No data for this element
c	photo	r	r: Cross-section in cm**2

c Tabulated fit coefficients - coeff(i,j,k) is the ith coefficient for the
c jth energy range for the kth element. i goes from 1 to 9 and j from 1 to
c 4.

      INTEGER MAXEDG, MAXCOF, MAXZ
      PARAMETER(MAXEDG=3, MAXCOF=9, MAXZ=79)

      DOUBLE PRECISION coeff(MAXCOF, MAXEDG+1, MAXZ)
      DOUBLE PRECISION edge(MAXEDG, MAXZ)
      DOUBLE PRECISION Elog, E1, E2, X, Xt

      REAL f(MAXZ)

      INTEGER ncoefs(MAXEDG+1, MAXZ)
      INTEGER nedges(MAXZ)
      INTEGER i, j, i1, i2, k

      INTEGER i3, i4

      REAL helxsc, helyan
      EXTERNAL helxsc, helyan

c ***** Coefficients for polynomial fits to cross-sections *****

c Hydrogen

      DATA ((coeff(i,j,1),i=1,9),j=1,4)/
     &  21.46941, 0.9398479, -0.1492932, 5.4634294d-3,
     &  32*0./

c Helium

      DATA ((coeff(i,j,2),i=1,9),j=1,4)/
     &  14.61546, 4.682793, -0.7323856, 4.6526663d-2, -1.1172282d-3,
     &  31*0./

c Lithium

      DATA ((coeff(i,j,3),i=1,9),j=1,4)/ 36*0. /

c Beryllium

      DATA ((coeff(i,j,4),i=1,9),j=1,4)/ 36*0. /

c Boron

      DATA ((coeff(i,j,5),i=1,9),j=1,4)/ 36*0. /

c Carbon

      DATA ((coeff(i,j,6),i=1,9),j=1,4)/
     &  8.74161, 7.13348, -1.14604, 0.0677044, 5*0.,
     &  3.81334, 8.93626, -1.06905, 0.0422195,
     &  23*0./

c Nitrogen

      DATA ((coeff(i,j,7),i=1,9),j=1,4)/
     &  9.24058, 7.02985, -1.08849, 0.0611007, 5*0.,
     &  -13.0353, 15.4851, -1.89502, 0.0769412,
     &  23*0./

c Oxygen

      DATA ((coeff(i,j,8),i=1,9),j=1,4)/
     &  2.57264, 10.9321, -1.79383, 0.102619, 5*0.,
     &  16.53869, 3.6428144, -0.3177744, 7.9471897d-3,
     &  23*0./

c Fluorine

      DATA ((coeff(i,j,9),i=1,9),j=1,4)/ 36*0. /

c Neon

      DATA ((coeff(i,j,10),i=1,9),j=1,4)/
     &  -3.04041, 13.0071, -1.93205, 0.0977639, 5*0.,
     &  17.6007, 3.29278, -0.263065, 5.68290d-3,
     &  23*0./

c Sodium

      DATA ((coeff(i,j,11),i=1,9),j=1,4)/
     &  -2737.598, 2801.704, -1009.892, 87.16455, 43.20644,
     &  -15.27259, 2.180531, -0.1526546, 4.3137977d-3,
     &  1.534019, 9.261744, -0.9914126, 3.5278253d-2,
     &  23*0./

c Magnesium

      DATA ((coeff(i,j,12),i=1,9),j=1,4)/
     &  7.107172, 3.7359418, 7*0.,
     &  -81.32915, 65.2775, -15.00826, 1.558686, -6.1339621d-2, 4*0.,
     &  -9.161526, 13.07448, -1.435878, 5.2728362d-2,
     &  14*0./

c Aluminium

      DATA ((coeff(i,j,13),i=1,9),j=1,4)/
     &  26.90487, -6.135221, 1.175546, 6*0.,
     &  -38.1232, 29.5161, -4.45416, 0.226204, 5*0.,
     &  14.6897, 4.22743, -0.344185, 8.18542d-3,
     &  14*0./

c Silicon

      DATA ((coeff(i,j,14),i=1,9),j=1,4)/
     &  -3.066295, 10.006248, -0.9627411, 6*0.,
     &  -182.7217, 128.061, -29.47269, 3.03284, -0.1173096, 4*0.,
     &  -33.39074, 21.42992, -2.385117, 8.887583d-2,
     &  14*0./

c Phosphorus

      DATA ((coeff(i,j,15),i=1,9),j=1,4)/ 36*0. /

c Sulphur

      DATA ((coeff(i,j,16),i=1,9),j=1,4)/
     &  598.2911, -675.2265, 308.1133, -68.99324, 7.62458,
     &  -0.3335031, 3*0,
     &  3994.831, -3690.886, 1417.287, -287.9909, 32.70061,
     &  -1.968987, 4.9149349d-2, 2*0.,
     &  -22.49628, 17.24599, -1.848444, 6.6506132d-2,
     &  14*0./

c Chlorine

      DATA ((coeff(i,j,17),i=1,9),j=1,4)/
     &  6253.247, -8222.248, 4491.675, -1302.145, 211.4881,
     &  -18.25547, 0.6545154, 2*0.,
     &  -233.0502, 146.9776, -31.12463, 2.938618, -0.104096,
     &  4*0.,
     &  -23.74675, 17.50997, -1.857953, 6.6208832d-2,
     &  14*0./

c Argon

      DATA ((coeff(i,j,18),i=1,9),j=1,4)/
     &  -330.3509, 270.7433, -78.90498, 10.35983, -0.5140201,
     &  4*0.,
     &  -5.71870, 8.85812, -0.307357, 0.00169351, -0.0138134,
     &  0.00120451, 3*0,
     &  19.1905, 2.74276, -0.164603, 0.00165895,
     &  14*0./

c Potassium

      DATA ((coeff(i,j,19),i=1,9),j=1,4)/ 36*0. /

c Calcium

      DATA ((coeff(i,j,20),i=1,9),j=1,4)/
     &  -873.972, 868.5231, -339.678, 66.83369, -6.590398,
     &  0.2601044, 3*0,
     &  -3449.707, 2436.409, -682.0668, 95.3563, -6.655018,
     &  0.1854492, 3*0,
     &  18.89376, 2.709646, -0.1377201,
     &  15*0./

c Scandium

      DATA ((coeff(i,j,21),i=1,9),j=1,4)/ 36*0. /

c Titanium

      DATA ((coeff(i,j,22),i=1,9),j=1,4)/ 36*0. /

c Vanadium

      DATA ((coeff(i,j,23),i=1,9),j=1,4)/ 36*0. /

c Chromium

      DATA ((coeff(i,j,24),i=1,9),j=1,4)/
     &  -0.4919405, 15.66939, -5.199775, 1.086566, -0.1196001,
     &  5.2152011d-3, 3*0,
     &  27.29282, 0.2966640, 7*0.,
     &  -15.2525, 16.23729, -1.966778, 8.062207d-2, 5*0.,
     &  8.307041, 5.008987, -0.2580816, 6*0./

c Manganese

      DATA ((coeff(i,j,25),i=1,9),j=1,4)/ 36*0. /

c Iron

      DATA ((coeff(i,j,26),i=1,9),j=1,4)/
     &  -15.07332, 21.94335, -4.862457, 0.5573765, -3.0065542d-2,
     &  4.9834867d-4, 3*0,
     &  -253.0979, 138.4238, -25.47119, 2.08867, -6.4264648d-2,
     &  4*0.,
     &  -1.037655 , 7.022304, -0.3638919,
     &  15*0./

c Cobalt

      DATA ((coeff(i,j,27),i=1,9),j=1,4)/
     &  9.171919, 3.5721176, 7*0.,
     &  -6.910097, 13.58385, -1.873453, 9.1612935d-2, 5*0.,
     &  13.96877, 2.128918, 0.1149042, 4.9106661d-02, -1.4725224d-02,
     &  8.3086651d-04, 3*0.,
     &  28.72910, 0.4456830, 7*0./

c Nickel

      DATA ((coeff(i,j,28),i=1,9),j=1,4)/
     &  -7.919931, 14.06475, -1.935318, 9.3929626d-2, 5*0.,
     &  3.71129, 8.45098, -0.896656, 0.0324889, 5*0.,
     &  28.4989, 0.485797,
     &  16*0./

c Gold   !! Important : this is only valid above 220 eV !!

      DATA ((coeff(i,j,79),i=1,9),j=1,4)/
     &  -27.40668, 18.11780, -1.869548, 6.2878355D-02, 5*0.,
     &  26.50165, 0.6361220, 7*0.,
     &  -33.83069, 19.80218, -1.989242, 6.7341216D-02,
     &  14*0./

c Number of coefficients (-1 => range does not exist)
c H - Ne

      DATA ((ncoefs(i,j),i=1,4),j=1,10)/
     &  4, -1, -1, -1,
     &  5, -1, -1, -1,
     &  -1, -1, -1, -1,
     &  -1, -1, -1, -1,
     &  -1, -1, -1, -1,
     &  4, 4, -1, -1,
     &  4, 4, -1, -1,
     &  4, 4, -1, -1,
     &  -1, -1, -1, -1,
     &  4, 4, -1, -1/

c Na-Ar

      DATA ((ncoefs(i,j),i=1,4),j=11,18)/
     &  9, 4, -1, -1,
     &  2, 5, 4, -1,
     &  3, 4, 4, -1,
     &  3, 5, 4, -1,
     &  -1, -1, -1, -1,
     &  6, 7, 4, -1,
     &  7, 5, 5, -1,
     &  5, 6, 4, -1/

c K-Ni

      DATA ((ncoefs(i,j),i=1,4),j=19,28)/
     &  -1, -1, -1, -1,
     &  6, 6, 3, -1,
     &  -1, -1, -1, -1,
     &  -1, -1, -1, -1,
     &  -1, -1, -1, -1,
     &  6, 2, 4, 3,
     &  -1, -1, -1, -1,
     &  6, 5, 3, -1,
     &  2, 4, 6, 2,
     &  4, 4, 2, -1/

c Cu-Pt

      DATA ((ncoefs(i,j),i=1,4),j=29,78) / 200*-1/

c Au

      DATA (ncoefs(i,79),i=1,4) / 4, 2, 4, -1/


c *****  Edge energies (in eV)  *****

c H-Ne

      DATA ((edge(i,j),i=1,3),j=1,10) /
     & 3*1.d32,
     & 3*1.d32,
     & 3*1.d32,
     & 3*1.d32,
     & 3*1.d32,
     & 284.0, 2*1.d32,
     & 401.0, 2*1.d32,
     & 531.7, 2*1.d32,
     & 3*1.d32,
     & 867.0, 2*1.d32/

c Na-Ar

      DATA ((edge(i,j),i=1,3),j=11,18) /
     & 1071.7, 2*1.d32,
     & 49.45, 1303.4, 1.d32,
     & 72.78, 1559.9, 1.d32,
     & 100.6, 1840.0, 1.d32,
     & 3*1.d32,
     & 165.0, 2470.5, 1.d32,
     & 202.0, 2819.6, 1.d32,
     & 245.0, 3202.9, 1.d32/

c K-Ni

      DATA ((edge(i,j),i=1,3),j=19,28) /
     &  3*1.d32,
     &  349.31, 4038.1, 1.d32,
     &  3*1.d32,
     &  3*1.d32,
     &  3*1.d32,
     &  598.0, 691.0, 5988.8,
     &  3*1.d32,
     &  707.4, 7111.2, 1.d32,
     &  61., 793.8, 7709.5,
     &  853.6, 8331.6, 1.d32/

c Cu-Pt

      DATA ((edge(i,j),i=1,3),j=29,78) / 150*1.d32/

c Au

      DATA (edge(i,79),i=1,3) / 2220., 2743.9, 1.d32/

c Number of edges (-1 => no tabulated data)

      DATA nedges /
     &  0, 0, -1, -1, -1, 1, 1, 1, -1, 1,
     &  1, 2, 2, 2, -1, 2, 2, 2,
     &  -1, 2, -1, -1, -1, 3, -1, 2, 3, 2,
     &  50*-1, 2/

c cm*cm/g to barns conversion factor

      DATA f/
     *  1.674, 6.646, 0, 0, 0, 19.94, 23.26, 26.56, 0, 33.50,
     *  38.17, 40.35, 44.80, 46.63, 0, 53.24, 58.86, 66.33, 0, 66.54,
     *  0, 0, 0, 86.33, 0, 92.72, 97.85, 97.48, 50*0., 327.0/

      photo = 0.

c Special case for Helium

      IF (Z .EQ. 2) THEN
         IF ( versn .EQ. 2 ) THEN
            photo = helxsc(500*(keV1+keV2)) * EXP(-55.26204) * f(Z)
         ELSEIF ( versn .EQ. 3 ) THEN
            photo = helyan(500*(keV1+keV2)) * EXP(-55.26204) * f(Z)
         ENDIF
         RETURN
      ENDIF

c Other elements

      IF (nedges(Z) .EQ. -1) THEN
         status = 1
         photo = 0.
         RETURN
      ELSE
         status = 0
      ENDIF

      E1 = 1000*keV1

      E2 = 1000*keV2

c Find the appropriate range to contain the lower energy.

      i1 = 1
c      DO WHILE ( (E1.GE.edge(i1,Z)) .AND. (i1.LE.nedges(Z)) )
c         i1 = i1 + 1
c      ENDDO

	DO i3=1,3
		IF ((i1.LE.nedges(Z))) THEN
                   IF (E1.GE.edge(i1,Z)) THEN 
			i1 = i1 + 1
			ENDIF
                   ENDIF
		ENDDO


c Find the appropriate range to contain the upper energy.

c      i2 = i1
c      IF (i2.LE.nedges(Z)) THEN
c         DO WHILE ( (E2.GE.edge(i2,Z)) .AND. (i2.LE.nedges(Z)) )
c            i2 = i2 + 1
c         ENDDO
c         ENDIF

      i2 = i1 
	DO i4=1,3
		IF ((i2.LE.nedges(Z))) THEN
                   IF (E2.GE.edge(i2,Z)) THEN 
			i2 = i2 + 1
			ENDIF
                   ENDIF
		ENDDO


c If these are the same then just sum up the cross-section
c for the midpoint of the bin.

      IF (i1 .EQ. i2) THEN

         Elog = log((E1+E2)/2)

         X = coeff(ncoefs(i1,Z),i1,Z)
         DO j = ncoefs(i1,Z)-1, 1, -1
            X = X * Elog + coeff(j,i1,Z)
         ENDDO

      ELSE

c First do the lower energy up to the first edge in the bin

         Elog = log((E1+edge(i1,Z))/2)

         Xt = coeff(ncoefs(i1,Z),i1,Z)
         DO j = ncoefs(i1,Z)-1, 1, -1
            Xt = Xt * Elog + coeff(j,i1,Z)
         ENDDO
         X = Xt * (edge(i1,Z)-E1)/(E2-E1)

c Now calculate the last edge in the bin up to the upper energy

         Elog = log((E2+edge(i2-1,Z))/2)

         Xt = coeff(ncoefs(i2,Z),i2,Z)
         DO j = ncoefs(i2,Z)-1, 1, -1
            Xt = Xt * Elog + coeff(j,i2,Z)
         ENDDO
         X = X + Xt * (E2-edge(i2-1,Z))/(E2-E1)

c Now add in any bits between edges in the bin

         DO k = i1+1, i2-1

            Elog = log((edge(k,Z)+edge(k-1,Z))/2)
            Xt = coeff(ncoefs(k,Z),k,Z)
            DO j = ncoefs(k,Z), 1, -1
               Xt = Xt * Elog + coeff(j,k,Z)
            ENDDO
            X = X + Xt * (edge(k,Z)-edge(k-1,Z))/(E2-E1)

         ENDDO

      ENDIF

c Do the exponential, put in the E**3 factor, and convert to cm**2.

      Elog = log((E1+E2)/2)
      photo = EXP(X-3*Elog-55.26204)*f(Z)

      RETURN
      END

c---------------------------------------------------------------------
      FUNCTION HELXSC(E)

      REAL HELXSC, E

C
C     Real Funcion : HELIUM
C     Source : Marr, G. V., and West, J. B., Atomic and Nuclear Data Tables,
C                (1976) 18, 497.
C             Oza, D. H., (1986), Phys. Rev. A, 33,  824.
C             Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
C                J. Phys. B., 20, 6457.
C
C     Description :
C     calculates mass absorption coefficient (mu/rho) in cm2/g for neutral
C     helium for the given energy in eV.
C     Cross sections come from experimental data compiled by Marr and
C     West (Atomic Data and Nuclear Data Tables (1976) 18, 497).
C     The four strongest autoionization resonances are taken into account;
C     numbers come from Oza (Phys Rev A (1986), 33, 824), and Fernley et al.
C     (J. Phys B (1987) 20, 6457).
C
C     Deficiencies :
C     works in the energy range from 30 eV to 10,000 eV

C     Bugs :
C     if any are found please report to the authors
C
C     History :
C     this subroutine replaces the previous version of HELIUM which
C     calculated mass absoprtion coefficients based on Henke's data
C     (Henke, B. L., et al., (1982), Atomic and Nuclear Data Tables, 27, 1).
C     This version of HELIUM returns mass  absorption coefficients which
C     are in better agreement with the best experiments as well as
C     theoretical models (see Chen, W. F., Cooper, G., and Brion, C. E.,
C     (1991), Phys. Rev. A, 44, 186).  This fortran-77 version of the
C     subroutine is based on Pat Jelinsky's program written in C
C     (obtained from EUVE Archive)
C
C     History :
C     04 jan 93 : original (19775::MBC)
C
C     23 feb 93 : comments added and modified to remove VAX
C                    fortran 77 extensions
C
C     21 sep 93 : further remnants of VAX fortran 77 extensions
C                    have been removed (19775::MBC)
C
C     23 sep 93 : bug in the FANO routine has been removed (19775::MBC)
C
C     Usage : FUNCTION HELIUM(E)
C            E = Energy in eV
C
C     Common Blocks :
C           none
C
C     Implicit :
C           none
C
C     Functions called by HELIUM
C           FANO
C
C------------------------------------------------------------------------------

C  Avogadro's number

      REAL AV
      PARAMETER (AV=6.022045E23)

C  atomic weight of hydrogen

      REAL AW
      PARAMETER (AW=4.0026E0)

      INTEGER IP, IF
      PARAMETER (IP=8, IF=4)

      REAL C1(IP), C2(IP), Q(IF), NU(IF), GAMMA(IF)
      REAL LAMBDA, X, Y, SIGMA, EPS

      INTEGER I

C polynomial coefficients for Marr and West data
      DATA C1 /-2.953607E1, 7.083061E0, 8.678646E-1, -1.221932E0,
     +         4.052997E-2, 1.317109E-1, -3.265795E-2, 2.500933E-3/

C polynomial coefficients for Marr and West data )
      DATA C2 /-2.465188E1, 4.354679E0, -3.553024E0, 5.573040E0,
     +          -5.872938E0, 3.720797E0, -1.226919E0, 1.576657E-1/

C parameters Q for resonances (Fernley et al. 1987)
      DATA Q /2.81E0, 2.51E0, 2.45E0, 2.44E0/

C parameters NU for resonances (Oza 1986)
      DATA NU /1.610E0, 2.795E0, 3.817E0, 4.824E0/

C parameters GAMMA for resonances (Oza 1986)
      DATA GAMMA /2.64061E-3, 6.20116E-4, 2.56061E-4,
     +             1.320159E-4/


C Calculate wavelength

      LAMBDA=12398.54E0/E
      X=ALOG10(LAMBDA)

c If > 503.97 then no absorption

      IF(LAMBDA.GT.503.97E0) THEN
         HELXSC=0.E0
         RETURN

c If < 46 then use first polynomial fit

      ELSEIF(LAMBDA.LT.46.E0) THEN
         Y=0.E0
         DO I=1,IP
            Y=Y + C2(I)*(X**(I-1))
         ENDDO

c Otherwise use second polynomial fit and include autoionization
c resonances

      ELSE

         Y=0.E0
         DO I=1,IP
            Y=Y + C1(I)*(X**(I-1))
         ENDDO

         EPS=911.2671E0/LAMBDA
         DO I=1,IF
            X=2.0*(EPS - 3.0E0 + 1.E0/(NU(I)*NU(I)) - 1.807317)/GAMMA(I)
            Y=Y + ALOG10((X-Q(I))*(X-Q(I))/(1.0E0 + X*X))
         ENDDO

      ENDIF

      SIGMA=10.E0**Y
      HELXSC=SIGMA*AV/AW
      END

C-----------------------------------------------------------------------------
C
C     Real Funcion : HELYAN
C     Source : Yan et al 1998, ApJ 496, 1044.
C             Oza, D. H., (1986), Phys. Rev. A, 33,  824.
C             Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
C                J. Phys. B., 20, 6457.
C
C     Description :
C     calculates mass absorption coefficient (mu/rho) in cm2/g for neutral
C     helium for the given energy in eV.
C     Cross sections come from a theoretical adjustment of experimental data
C     (from Samson et al 1994a and others) by Yan et al. (ApJ 1997).
C
C     The four strongest autoionization resonances are taken into account;
C     numbers come from Oza (Phys Rev A (1986), 33, 824), and Fernley et al.
C     (J. Phys B (1987) 20, 6457).
C
C     Deficiencies :
C     works in the energy range from 30 eV to 10,000 eV

C     Bugs :
C     if any are found please report to the authors
C
C     History :
C     1991: calculated mass absoprtion coefficients based on Henke's data
C     (Henke, B. L., et al., (1982), Atomic and Nuclear Data Tables, 27, 1).
C     1993: modified to include autoionization resonances.  The fortran-77
C     version of the subroutine was based on Pat Jelinsky's C program.
C     (obtained from EUVE Archive)
C     1997: realized Marr and West (1976:
C     cross sections are a poor match to current best estimates above 80 eV.
C     Changed continuum portion to cross sections recommended by Yan et al
C     (ApJ 1977)			Autoionization resonances were
C     retained as implemented by Jelinsky in the EUVE archive.
C
C
C     Usage : FUNCTION HELYAN(E)
C            E = Energy in eV
C
C     Common Blocks :
C           none
C
C     Implicit :
C           none
C
C     Functions called by HELYAN
C           FANO
C
C------------------------------------------------------------------------------
       FUNCTION HELYAN(E)

C    Type definitions :
C     IMPLICIT NONE
C    Global variables :
C    Structure definitions :
C    Function declarations :
       REAL FANO
C    Local constants :
       INTEGER IP
C         ( index through loop )
       PARAMETER (IP=6)
       INTEGER IF
C         ( index through loop )
       PARAMETER (IF=4)
       REAL AV
C         ( Avogadro's number )
       PARAMETER (AV=6.022045E23)
       REAL AW
C         ( atomic weight of hydrogen )
       PARAMETER (AW=4.0026E0)
C    Local variables :
        REAL LAMBDA
C          ( wavelength in Angstroms)
        REAL X
        REAL Y
        REAL SIGMA
C          ( cross section in cm2/atom)
        INTEGER I
C          ( index trough loop)
C     Import :
        REAL E
C          ( energy in eV)
C     Export :
        REAL HELYAN
C          ( cross section in cm**2/g)
C    Local data :
        REAL C1(IP)
        REAL EION
        REAL Q(IF)
        REAL NU(IF)
        REAL GAMMA(IF)

C          ( polynomial coefficients for Yan et al data)
        DATA C1 /-4.7416, 14.8200, -30.8678, 37.3584, -23.4585, 5.9133/

C          ( ionization edge in eV:)
        DATA EION  /24.58/

C          ( parameters Q for resonances (Fernley et al. 1987) )
        DATA Q /2.81E0, 2.51E0, 2.45E0, 2.44E0/

C          ( parameters NU for resonances (Oza 1986) )
        DATA NU /1.610E0, 2.795E0, 3.817E0, 4.824E0/

C          ( parameters GAMMA for resonances (Oza 1986) )
        DATA GAMMA /2.64061E-3, 6.20116E-4, 2.56061E-4,
     +             1.320159E-4/


C     Start :

        LAMBDA=12398.54E0/E

        X = E/EION

        IF(LAMBDA.GT.503.97E0) THEN
          HELYAN=0.E0
        ELSE
          Y=1.E0
          DO 2 I=1,IP
  2         Y=Y + C1(I)/(X**(I/2.))

C        Yan et al. cross section in cm**2/atom:
          SIGMA = 733.E-24/(1.E-3*E)**3.5 * Y

C        Add in autoionization resonances:
          DO 3 I=1,IF
  3          SIGMA = SIGMA * FANO(Q(I),NU(I),GAMMA(I),LAMBDA)

          HELYAN=SIGMA*AV/AW

        ENDIF

        RETURN
        END

C----------------------------------------------------------------------------
C
        FUNCTION FANO(A,B,C,LAMBDA)

C    Type definitions :
C     IMPLICIT NONE
C    Global variables :
C    Structure definitions :
C    Function declarations :
C    Local constants :
C    Local variables :
        REAL EPS
C          ( energy in Rydbergs )
        REAL EPSI
        REAL X
C          ( log_10 of wavelength in Angstroms )
C     Import :
        REAL A
C          ( Q coefficient (Fernley et al. 1987) )
        REAL B
C          ( NU coefficient (Oza 1986) )
        REAL C
C          ( GAMMA coefficient (Oza 1986) )
        REAL LAMBDA
C          ( wavelength in Angstroms )
C     Export :
        REAL FANO
C    Start :

        EPS=911.2671E0/LAMBDA
        EPSI=3.0E0 - 1.E0/(B*B) + 1.807317
        X=2.0*(EPS - EPSI)/C
        FANO=(X-A)*(X-A)/(1.0E0 + X*X)

        END
C-----------------------------------------------------------------------
