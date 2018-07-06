C     Least square fit of charge distribution parameters to elastic
C     electron scattering cross sections
C
C     Version 4C
C     This version was created in December of 1981 to run on
C     the Saclay NORD-500, a 32-bit machine.  The character handling
C     has been modified to conform to Fortran-77 standards.  The
C     program has been run through the CERN library program TIDY so
C     all statement numbers have changed.  All floating point variables
C     constants, and functions have been changed to double precision.
C     The new calls to the charge density routines and the differential
C     equation integration routine have been put in. Despite all of
C     these changes, the version number is still 4C, the version of
C     the current CYBER version, because nothing big has been altered.
C
C     Version 5A
C     In version 5A, an attempt has been made to use a structured
C     programming style to produce a more readable, changable, and
C     bug-free program (using NORD Fortran and Fortran 77 constructs).
C     The ability to fit generalized radial moments as well as cross
C     sections has been added. Physical constants have been updated.
C     In order to reduce overhead, routines have been added to allow
C     the storage of fitted data in (virtual) memory (CONDITIONAL-
C     COMPILATION V) rather than in temporary disk files (CONDITIONAL-
C     COMPILATION D). The definition of percent deviation has been
C     changed to conform to the Saclay standard. The ability to output
C     deviations from the fitted result in units of error bars has
C     been added.  Finally, large parts of the listing output can
C     now be supressed.
C
C     In version 5B, the program has been altered to conform as close
C     as possible to the Fortran 77 conventions. Special NORD Fortran
C     constructs such as DOFOR, ENDDO, and the -1:-1 in strings have
C     been avoided.  A bug in the plots has been fixed.  This version
C     runs on the Illinois Perkin-Elmer 3220 under OS32 (MTM).
C
C     Version 5B was converted to run on a DEC VAX running VMS.
C     Most alterations involved minor changes in I/O (like unit
C     numbers).  Note that this must compiled and run using the
C     G_FLOATING option. This option has been included explicitly
C     in the source.  The program was modified to work with
C     negative Z's for positron scattering.
C
C     In version 5C, temporary files are now scratch files so they
C     are deleted at the end of a run.  Also, the Coulomb functions
C     file may be specified. If the file exists, then Coulomb functions
C     are read from it.  If the file does not exist, it is created and
C     calculated cross sections are written to it.  If no name is
C     entered, then Coulomb functions are calculated but not stored.
C     Some improvements to the comments below have been made and a
C     bug in the cpu/clock time printout was fixed.
C
C     Version 5D includes the ravpro.inc parameter file in order to set
C     the default maximum number of phase shifts (NFMX).  This change
C     is required if version 1D or greater of RAVPRO is used.
C
C     Version 5E has been updated to use some standard (wlib) routines
C     and cpp constructs for conditional compilation.
C
C     Input data:
C     All input data is entered in free format except the first card.
C
C     Card 1: TITLE
C          TITLE:  A title for the fit which appears on printout and
C                  plots.
C
C     Card 2: Z0,A,MODEL,NITMAX,INRM,NGENM,CONV,STEP,NFL
C          Z0:     If not equal to 0, Z0 is the Z for the nucleus to
C                  be fit.  In this case, the first parameter of the
C                  charge density is taken to be an overall density
C                  normalization.  Its initial value is calculated
C                  to force the integral of the density to yield Z0.
C                  If equal to 0, the charge density routine is assumed
C                  to be correctly normalized.  Its integral, we assume,
C                  gives the Z of the nucleus and the identification
C                  of the first parameter as a normalization is not
C                  made.
C          A:      A for the nucleus to be fit.
C          MODEL:  The model number for the charge density for the
C                  nucleus to be fit. For linear models, like Fourier-
C                  Bessel and Lageure, this is the number of terms in
C                  the series to consider.
C          NITMAX: The maximum number of charge density fitting
C                  iterations to do. During such iterations, only the
C                  parameters of the charge density are varied.
C                  If NITMAX is 0 and INOFIT is also 0, then errors and
C                  CHI-square will be calculated with the initial
C                  parameter values but no fitting will be done.
C          INRM:   The number of normalization fitting iterations to
C                  do.  During such iterations, only the cross section
C                  normalizations are varied.  NITMAX charge density
C                  iterations are done before each normalization fitting
C                  iteration.
C          NGENM:  The number of generalized moments to be fit in
C                  addition to the cross sections.  This number must
C                  be less than the value of the parameter NGENMX which
C                  is 20 at the moment.
C          CONV:   Minimize size of phase shift to be considered during
C                  calculation of cross sections by RAVPRO (default is
C                  1E-6.
C          STEP:   Step size to use when integrating the Dirac equation
C                  inside the "fitting-on" radius.  (default is .05
C                  in units of energy/h-barc???)
C          NFL:    Number of terms in the legendre series to consider.
C                  (default is 50)
C
C     Card 3: IOLD,IABSP,IDEVP,IRMS,IZ,INOFIT,IRHOP,IDERB,IFFP,ICHP,ILPR
C          IOLD    0 => Cross sections are relative.  A set of cross
C                  sections for a reference nucleus is to be read in
C                  and, using these values and a model for the reference
C                  nucleus, the cross sections for the nucleus to be
C                  fit are to be converted to absolute cross sections.
C                  non-0 => Cross sections are absolute.  No reference
C                  nucleus data is to be read in.
C          IABSP   0 => Do not write absolute cross sections to the
C                  plot/punch output file.
C                  non-0 => Absolute cross sections are to be written to
C                  the plot/punch output file.
C          IDEVP   0 => No percent deviations of the data from the fit
C                  are to be written to the plot/punch output file.
C                  1 => Percent deviations of the data from the fit are
C                  to be written to the plot/punch output file.
C          IRMS    0 => No Lagrange multiplier constraint on the RMS
C                  radius.
C                  1 => RMS is constrained during the fit to remain at
C                  its initial value.
C          IZ      0 => No Lagrange multiplier constraint on the Z of
C                  the nucleus during the fit.  (generally this is
C                  not a very useful option since Z is usually known
C                  very accurately)
C                  1 => Z is constrained to be its initial value
C                  during the fit.
C          INOFIT: 0 => Fitting of the data is to be done.
C                  non-0 => No fitting of the data.  Absolute cross
C                  sections may be calculated though.
C          IRHOP:  0 => Charge density, and charge density error
C                  should not be written to the plot/punch output file.
C                  non-0 => Charge density, and charge density error
C                  should be written to the plot/punch output file after
C                  the fit.
C          IDERB:  0 => During the fitting process, initially, 1-point
C                  derivitives of the cross sections with respect to
C                  the parameters should be calculated and automatically
C                  this should change to two point calculation when
C                  fitting starts to diverge (near a minimum).
C                  non-0 => Always use 2-point calculation of deriv-
C                  itives.
C          IFFP:   0 => No form factor plotting information should be
C                  written to the plot/punch output file.
C                  1 => Experimental form factors and fitted form
C                  factors between 1 and 179 degrees at maximum exp-
C                  erimental energy should be written to the plot/
C                  punch output file.
C                  not 0 or 1 => Experimental form factors and fitted
C                  form factors between qeff=0 and qeff=1.5 times the
C                  maximum experimental qeff should be written to the
C                  plot/punch output file.
C          ICHP    0 => No deviations of the data from the fit in units
C                  of error bars are to be written to the plot/punch
C                  output file.
C                  1 => Deviations of the data from the fit in units of
C                  error bars are to be written to the plot/punch output
C                  file.
C          ILPR    0 => Print only new parameters and afew other things
C                  after each iteration of the fit.
C                  1 => Print detailed point-by-point analysis error
C                  matrix, and charge density after every iteration.
C
C     Card 4: (NPRINT(I),I=1,10),(NPUNCH(I),I=1,5)
C          NPRINT(I): Causes the following debugging information to be
C                  printed when set to a non-0 value. If requested
C                  during a fit, this can generate alot of output so
C                  be careful.
C                  1.  Coulomb function calculation subroutine (couls)
C                      input
C                  2.  Ravenhall program (RAVPRO) input
C                  3.  Charge distribution information (rho, v, norm,..)
C                  4.  Dirac equation solutions
C                  5.  Phase shifts
C                  6.  energy, angles, sigma in c/m, form factors
C                  7.  COULS integration results
C                  8.  Intermediate results of integration and angle
C                      dependant calculation
C                  9.  Intermediate integration results
C                 10.  folding integration results
C          NPUNCH(I): Causes the follownig debugging information to be
C                  written to the plot/punch output file when set to a
C                  non-0 value. If requested during a fit, this can
C                  generate alot of output so be careful.
C                  1.  Coulomb phase shifts
C                  2.  Amplitudes
C                  3.  Cross sections and angles
C                  4.  Unused
C                  5.  Unused
C
C     A card for each charge density parameter: PAR(I),EPS(I),DPAR(I)
C          PAR(I): The value of the parameter. The order and number
C                  of parameters depends on the model.
C          EPS(I): If non-0, EPS(I) is taken to be a minimum percent
C                  change in the corresponding parameter between two
C                  iterations.  If all parameters change by less than
C                  the EPS values, fitting is stopped.
C                  If 0, the corresponding parameter is fixed during
C                  the fit to its initial value.
C          DPAR(I):The parameter step size to be used when calculating
C                  derivitives of the cross sections with respect to
C                  the parameter, entered as a percentage of the para-
C                  meter (default is .01%).  If the corresponding para-
C                  meter is 0, then the step size is set to .01.
C     Note that there will be a certain number of parameters required
C     by the charge density subroutine to define the density.  An
C     additional parameter is always added by FITR.  This parameter is
C     an overall cross section normalization (much like the group
C     normalizations discussed below).  This normalization multiplies
C     all theoretical cross sections.  It should not be fit (its values
C     of DPAR and EPS should be set to 0) if the group normalizations
C     are fit.
C
C     A card for each of the experimentally measured generalized
C     moments as defined by Ford and Rinker (Phys Rev C7 (1973), 1208)
C     There should be NGENM cards of this type. Generalized Moment
C     Card: GMOMEX,GMOMAL,RMOME,DRMOME,LEVNAM
C          GMOMEX: Exponent of r
C          GMOMAL: The constant ALPHA in the exponential
C          RMOME:  The experimental value of the moment (fm)
C          DRMOME: The absolute uncertainty in the experimental
C                  moment (fm)
C          LEVNAM: A twelve-character identifier for the radial moment
C                  (generally the atomic level from which it was
C                  derived)
C
C     What follows next are a number of blocks of experimental data
C     (at the moment a maximum of 50 are allowed).  For each block,
C     there is one "energy card" and up to 200 (or 66 if all cross
C     sections are folded over acceptance and/or multiple scattering
C     angle) "angle cards".
C     Energy Card: EX,NANG,NF,X0,OTEST,GNRM,ERCOM,ID
C          EX:     Incident beam energy.
C          NANG:   Number of angle cards to follow.
C          NF:     Maximum number of phase shifts (default???)
C          X0:     Matching on radius??? (default???)
C          OTEST:  Maximum angle for which small angle calculation
C                  will be done (in degrees). (default???)
C          GNRM:   Normalization for this block of cross sections.
C                  If 0, the normalization is taken to be the same as
C                  the previous block or 1 if there was no previous
C                  block.  If negative, the normalization is the
C                  absolute value of GNRM and it is fixed during
C                  normalization fitting iterations.
C          ERCOM:  A common (reproducability) error which is added
C                  in quadrature to all cross section errors, entered
C                  as a percentage of the cross sections.
C          ID:     An identifier of up to ten characters which is used
C                  in plots and print-out to identify the data set.
C
C     Angle Cards: THX,SX,SXE,DELTH,DELPH,XMS
C          THX:    Angle in degrees
C          SX:     Cross section in fm**2 per steradian
C          SXE:    Absolute Cross section error in fm**2 per steradian
C          DELTH:  Horizontal (in scattering plane) acceptance of the
C                  spectrometer for folding of calculated cross section
C                  over solid angle. 0 => no folding.  This is the half
C                  angle in degrees.
C          DELPH:  Vertical acceptance of the spectrometer for folding
C                  of the cross sections over solid angle. 0 => no
C                  folding. This is the half angle in degrees.
C          XMS:    If positive, XMS is taken to be the multiple scatter-
C                  ing angle in degrees for folding.  If negative, it
C                  is taken to be the target thickness in g/cm**2 and
C                  the multiple scattering angle is calculated.  0 => no
C                  folding over multiple scattering angle is to be done.
C
C     If IOLD was 0, information must be supplied for the reference
C     nucleus.  This consists if some data describing the model followed
C     by blocks of data measured at about the same energy and angles as
C     the data for the nucleus to be fit.
C     First Card: ZREF,AREF,MODELN
C          ZREF:   Z of the reference nucleus.
C          AREF:   A of the reference nucleus.
C          MODELN: Model number for the reference nucleus model.  For
C                  linear models, this is generally the number of
C                  terms in the series to consider.
C
C     A card for each charge density parameter PARAM
C          PARAM:  The value of a charge density parameter. The order
C                  and number of parameters depends on the model.
C
C     Reference nucleus energy cards: EXR,NFR,X0R,OTESTR
C          EXR:    The energy of the beam
C          NFR:
C          X0R:    The fitting-on radius at this energy (default is???)
C          OTESTR:
C
C     Reference nucleus angle cards (there must be NANG of these where
C     NANG is the value from the corresponding block of data for the
C     nucleus to be fit): THR,SR,RADIUS,DELTH,DELPH,XMS
C          THR:    The angle for the cross section in degrees.
C          SR:     The reference nucleus cross section.
C          RADIUS: The error in the reference nucleus cross section.
C                  (also used elsewhere in the program, thus the name)
C          DELTH:  Horizontal (in scattering plane) acceptance of the
C                  spectrometer for folding of calculated cross section
C                  over solid angle. 0 => no folding.
C          DELP:   Vertical acceptance of the spectrometer for folding
C                  of the cross sections over solid angle. 0 => no
C                  folding.
C          XMS:    If positive, XMS is taken to be the multiple scatter-
C                  ing angle in degrees for folding.  If negative, it
C                  is taken to be the target thickness in g/cm**2 and
C                  the multiple scattering angle is calculated.  0 => no
C                  folding over multiple scattering angle is to be done.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Include parameter file to define RAVPRO array sizes.
      INCLUDE 'ravpro.inc'
C     Version number of this program
      CHARACTER VER*(*)
      PARAMETER (VER='5E')
C     Last edit date for this program
      CHARACTER EDIT*(*)
      PARAMETER (EDIT='10/14/04')
C     The maximum number of parameters that can be fit
      PARAMETER (NPMX=29)
C     The dimensions of PARAM, PAR, DPAR, and the first dimension of
C     DERIV should be one greater than NPMX
      PARAMETER (NPMXP1=NPMX+1)
C     The dimensions of BETA, CURV, ALFA, ERROR, and EPS should be 3
C     larger than the number of density parameters, NPMX.  The second
C     to last element of these arrays corresponds to the target
C     thickness normalization constant; the next to last element is for
C     constraining the RMS radius; and the last element is for
C     constraining Z.
      PARAMETER (NPMXP3=NPMX+3)
C     The maximum number of normalization groups
      PARAMETER (NGRMX=50)
C     The maximum number of points in a given data set
      PARAMETER (NPTSMX=200)
C     The number of radius values used in integration and printout
C     of the charge density
      PARAMETER (NRAD=151)
C     The maximum of NPTSMX and NRAD
      PARAMETER (NRADMX=200)
C     The maximum number of generalized moments
      PARAMETER (NGENMX=20)
      CHARACTER TITLE*70,FILE*80,TI*8,DA*8,ID*10,IDS*10,RHONA*50
      CHARACTER RHONAR*50,LEVNAM(NGENMX)*12,OUTPUT*132,PRGNAM*80
      LOGICAL LAST,NEWIT,LAMDWN,DERIVB,ADERCM,COULFI,RDCOUL,LNGPRI
      INTEGER NPRINT(10),NPUNCH(5)
      REAL TIMCPU,TIMCLK
      DOUBLE PRECISION SR(NPTSMX),SX(NPTSMX),THX(NPTSMX),THR(NPTSMX)
      DOUBLE PRECISION SXE(NPTSMX),DELTH(NPTSMX),DELPH(NPTSMX)
      DOUBLE PRECISION ST(NPTSMX),XMS(NPTSMX),RADIUS(NRADMX),RHO(NRADMX)
      DOUBLE PRECISION RHOERR(NRADMX),FRN(NFMX),GRN(NFMX),FIN(NFMX)
      DOUBLE PRECISION GIN(NFMX),FRNR(NFMX),GRNR(NFMX),FINR(NFMX)
      DOUBLE PRECISION GINR(NFMX),SXT(NGRMX),STT(NGRMX),PARAM(NPMX)
      DOUBLE PRECISION PAR(NPMXP1),DPAR(NPMXP1),DERIV(NPMXP1,NPTSMX)
      DOUBLE PRECISION BETA(NPMXP3),CURV(NPMXP3,NPMXP3),EPS(NPMXP3)
      DOUBLE PRECISION ALFA(NPMXP3,NPMXP3),ERROR(NPMXP3,NPMXP3)
      DOUBLE PRECISION GMOMEX(NGENMX),GMOMAL(NGENMX),RMOME(NGENMX)
      DOUBLE PRECISION DRMOME(NGENMX)
      COMMON /CONST/ XMN,HBARC,ALPHA,XME,PI,RAD,RTD
      COMMON /COUL/ FRN,GRN,FIN,GIN
      COMMON /IO/ NPRINT,NPUNCH
      COMMON /DATA/ EX,NANG,NF,X0,OTEST,GNRM,GNORMS,GNORM,ERCOM,
     + ERCOMS,THX,SX,SXE,DELTH,DELPH,XMS,ID
      COMMON /CHARGE/ MODELN,PARAM
      DATA INULL/100916/,ILINE/100016/,IDOT/100316/,IEND/-1/
      DATA IPOINT/100900/
      DATA LAMDWN/.FALSE./,DERIVB/.FALSE./
C
      CALL TOD(TI,DA)
      CALL HELLOM(6,'FITR',VER,EDIT,PRGNAM)
      CALL PROMPT(6,'Data input file? ')
      READ (5,11) FILE
11    FORMAT (A)
#     ifdef linux
C         Linux g77 does not know about the "READONLY" attribute.
          OPEN (7,FILE=FILE,STATUS='OLD')
#     else
          OPEN (7,FILE=FILE,STATUS='OLD',READONLY)
#       endif
      CALL PROMPT(6,'Listing output file? ')
      READ (5,11) FILE
#     ifdef unix
          OPEN (8,FILE=FILE,STATUS='UNKNOWN')
#     else
C         We assume VMS if not Unix.
          OPEN (8,FILE=FILE,STATUS='NEW')
#       endif
      CALL PROMPT(6,'Plot/punch output file? ')
      READ (5,11) FILE
      LAST=FILE.NE.' '
      IF (LAST) THEN
#         ifdef unix
              OPEN (9,FILE=FILE,STATUS='UNKNOWN')
#         else
C             We assume VMS if not Unix.
              OPEN (9,FILE=FILE,STATUS='NEW')
#           endif
      ENDIF
#     ifdef linux
C         Linux g77 assumes scratch files have no names.
          OPEN (10,STATUS='SCRATCH',FORM='UNFORMATTED')
          OPEN (11,STATUS='SCRATCH',FORM='UNFORMATTED')
#     else
          CALL TMPNAM(FILE)
          OPEN (10,FILE=FILE,STATUS='SCRATCH',FORM='UNFORMATTED')
          CALL TMPNAM(FILE)
          OPEN (11,FILE=FILE,STATUS='SCRATCH',FORM='UNFORMATTED')
#     endif
      READ (7,30,END=1475) TITLE
30    FORMAT (A)
      WRITE (8,150) PRGNAM
150   FORMAT ('1',A)
      CALL NOSPAC(TITLE,I,J)
      WRITE (8,221) TITLE(1:J)
      READ (7,*,END=1475) Z0,A,MODEL,NITMAX,INRM,NGENM,CONV,STEP,NFL
      NGENM=MIN0(NGENM,NGENMX)
      IF (CONV.LE.0.D0) CONV=1D-6
      IF (STEP.LE.0.D0) STEP=5D-2
      IF (NFL.LE.0) NFL=50
      MODELN=MODEL
      RSQR=1.2D0*A**.33333333D0
      READ (7,*) IOLD,IABSP,IDEVP,IRMS,IZ,INOFIT,IRHOP,IDERB,IFFP,ICHP,
     + ILPR
      READ (7,*) (NPRINT(I),I=1,10),(NPUNCH(I),I=1,5)
      IF (.NOT.LAST) THEN
C         If no plot/punch file was specified, then supress all output
C         to the file
          IABSP=0
          IDEVP=0
          IRHOP=0
          IFFP=0
          ICHP=0
          DO 900 I=1,5
              NPUNCH(I)=0
900       CONTINUE
      ENDIF
      IF (INOFIT.EQ.0) THEN
          DERIVB=IDERB.NE.0
          CALL CHGDE3 (RHONA,NPARM1)
          NPAR=NPARM1+1
          NPARP1=NPAR+1
          NPARP2=NPARP1+1
          IF ((NPAR.LE.1).AND.(Z0.NE.0.D0))
     +     STOP 'Fitting not possible with this model'
          DO 901 I=1,NPAR
              READ (7,*) PAR(I),EPS(I),DPAR(I)
              PARAM(I)=PAR(I)
              IF (DPAR(I).EQ.0.D0) DPAR(I)=.1D0
              EPS(I)=EPS(I)*.01D0
              DPAR(I)=DPAR(I)*.01D0
901       CONTINUE
C         Determine if any of the model parameters are being fit
          LAST=.TRUE.
          IF (NPARM1.GE.1) THEN
              DO 902 I=1,NPARM1
                  LAST=LAST.AND.EPS(I).LE.0.D0
902           CONTINUE
          ENDIF
          IF (LAST) THEN
C             No model parameters fit so the RMS and Z can not be fixed
              IRMS=0
              IZ=0
C             If the over-all normalization is also not fit, then no fit
              IF (EPS(NPAR).LE.0.D0) NITMAX=0
          ENDIF
          EPS(NPARP1)=IRMS
          EPS(NPARP2)=IZ
          IF (PAR(NPAR).EQ.0.D0) PAR(NPAR)=1.D0
          IF (Z0.EQ.0.D0) THEN
C             Z (entered) equal to 0 means that model automatically
C             gives the value of Z
              Z=ZCALC(0)
          ELSE
C             Z (entered) not equal to 0 means first parameter of model
C             is an over all normalization of the density.
              PARAM(1)=1.D0
              Z=Z0
              PAR(1)=Z0/ZCALC(0)
          ENDIF
          WRITE (8,160) Z,A,RHONA,NITMAX,INRM
160       FORMAT ('0Elastic cross section fitting for Z = ',F6.1,
     +     ' , A = ',F6.1/'0The model is ',A/
     +     '0Fitting stopped after',I3,' iterations and',I3,1X,
     +     'normalizations or'/
     +     ' when the variation between iterations in all free',
     +     ' parameters'/' is less than the value listed below.'/
     +     '0Parameter',T20,'Initial',T35,'Min. variation to',T57,
     +     'Variation for'/T21,'value',T38,'continue fit',T59,
     +     'derivative'/1X,68('-'))
          DO 903 I=1,NPAR
C             We use output to a character string here to avoid multiple
C             lines in the output file which end up being printed on
C             one line.
              OUTPUT=' '
              WRITE (OUTPUT,170) I,PAR(I)
170           FORMAT (T6,I2,T15,1PE15.8)
              IF (EPS(I).NE.0.D0) THEN
                  WRITE (OUTPUT(39:45),190) EPS(I)
190               FORMAT (2PF6.4,'%')
                  WRITE (OUTPUT(58:64),190) DPAR(I)
              ENDIF
              CALL NOSPAC(OUTPUT,J,K)
              WRITE (8,221) OUTPUT(1:K)
903       CONTINUE
          IF (IRMS.NE.0) WRITE (8,200)
200       FORMAT ('0The RMS radius is constrained.')
          IF (IZ.NE.0) WRITE (8,210)
210       FORMAT ('0Z is constrained.')
          IF (Z0.NE.0.D0) WRITE (8,220)
220       FORMAT ('0The initial value of the charge density',
     +     ' normalization (parameter 1) will be calculated'/
     +     ' using the initial values of the other parameters')
          Z0=Z
          WRITE (8,230) NPAR,CONV,STEP,NFL
230       FORMAT ('0Target thickness normalization is parameter',I3/
     +     '0Calculation stopped for phase shifts of less than',
     +     1PE15.8/' Integration step size = ',1PE15.8/' NFL = ',I3)
      ELSE
          WRITE (8,250)
250       FORMAT (/'0No fitting done.')
      ENDIF
C
C     Read data for nucleus to be fit
      IF (NGENM.NE.0) THEN
          DO 904 I=1,NGENM
              READ (7,*) GMOMEX(I),GMOMAL(I),RMOME(I),DRMOME(I),
     +         LEVNAM(I)
904       CONTINUE
      ENDIF
      NENERG=0
      EXM=0.D0
      IF (IOLD.EQ.0) WRITE (8,270)
270   FORMAT ('1Experimental data'/' Relative cross sections for',
     + ' nucleus to be fit')
280   READ (7,*,END=310) EX,NANG,NF,X0,OTEST,GNRM,ERCOM,ID
      IF (EX.NE.0.0D0) THEN
          IF (X0.LE.0.0D0) X0=3.0D0*RSQR*EX/HBARC
          X0=STEP*INT(X0/STEP)
          IF (NF.LE.0) NF=1.0D0+3.0D0*X0
          NF=MIN(NF,NFMX-1)
          IF ((OTEST.LE.0.0D0).AND.(2.D0*EX*RSQR.LE.HBARC))
     +     OTEST=180.0D0
          IF (OTEST.LE.0.0D0) OTEST=2.0D0*RTD*DASIN(HBARC/
     +     (2.D0*EX*RSQR))
          IF (EX.GT.EXM) THEN
              EXM=EX
              X0M=X0
              NFM=NF
              OTESTM=OTEST
          ENDIF
          IF ((GNRM.EQ.0.0D0).AND.(NENERG.EQ.0)) GNRM=-1.0D0
          IF (GNRM.NE.0.D0) THEN
              GNORMS=GNRM
              ERCOMS=ERCOM*.01D0
              IDS=ID
              GNORM=DABS(GNORMS)
          ENDIF
          ERCOM=ERCOMS
          ID=IDS
          READ (7,*) (THX(J),SX(J),SXE(J),DELTH(J),DELPH(J),XMS(J),
     +     J=1,NANG)
          NENERG=NENERG+1
          IF (IOLD.EQ.0) CALL PDATA (NENERG,EX,NF,X0,OTEST,GNORM,0.D0,
     +     ERCOM,ID,THX,SX,SXE)
          CALL WDATA (NENERG)
          GO TO 280
      ENDIF
C
310   IF (IOLD.EQ.0) THEN
C         Change experimental relative cross sections to absolute
          READ (7,*) ZREF,AREF,MODELN
          CALL CHGDE3 (RHONAR,NPARR)
          IF (NPARR.GT.0) READ (7,*) (PARAM(I),I=1,NPARR)
          WRITE (8,320)
320       FORMAT ('1Absolute cross sections for reference nucleus')
          IF ((NPARR.LE.0).AND.(ZREF.NE.0.D0))
     +     STOP 'Z can not be determined for reference nucleus'
          IF (ZREF.LE.0.D0) THEN
              ZREF=ZCALC(0)
          ELSE
              PARAM(1)=1.D0
              PARAM(1)=ZREF/ZCALC(0)
          ENDIF
          WRITE (8,360) ZREF,AREF,RHONAR
360       FORMAT (/'0The reference nucleus is Z = ',F6.1,' , A = ',F6.1/
     +     ' The reference model is ',A)
          IF (NPARR.GT.0) WRITE (8,370) (PARAM(I),I=1,NPARR)
370       FORMAT (' with the following charge distribution parameters:'/
     +     ,5(1PE15.8))
          RSQRR=1.2D0*AREF**.333333333D0
          DO 905 I=1,NENERG
              CALL RDATA (I)
              READ (7,*) EXR,NFR,X0R,OTESTR
              IF (X0R.LE.0.0D0) X0R=3.0D0*RSQRR*EXR/HBARC
              X0R=STEP*INT(X0R/STEP)
              IF (NFR.LE.0) NFR=1.0D0+3.0D0*X0R
              NFR=MIN(NFR,NFMX-1)
              IF ((OTESTR.LE.0.0D0).AND.(2.D0*EXR*RSQRR.LE.HBARC))
     +         OTESTR=180.0D0
              IF (OTESTR.LE.0.0D0) OTESTR=2.0D0*RTD*DASIN(HBARC/
     +         (2.D0*EXR*RSQRR))
              READ (7,*) (THR(J),SR(J),RADIUS(J),DELTH(J),DELPH(J),
     +         XMS(J),J=1,NANG)
              CALL COULS (ZREF*ALPHA,X0R,NFR)
              CALL RAVPRO (ZREF*ALPHA,EXR,EXCM,AREF,CONV,X0R,STEP,
     +         OTESTR,THR,RHO,ST,RHOERR,NANG,NFR,NFL,DELTH,DELPH,
     +         XMS,0,RMS)
              DO 906 J=1,NANG
                  SXE(J)=DSQRT((SXE(J)/SX(J))**2+(RADIUS(J)/SR(J))**2)
                  SX(J)=SX(J)*ST(J)/SR(J)
                  SXE(J)=SXE(J)*SX(J)
906           CONTINUE
              CALL PDATA (I,EXR,NFR,X0R,OTESTR,0.D0,0.D0,0.D0,
     +         '          ',THR,SR,RADIUS)
              CALL WDATA (I)
905       CONTINUE
      ENDIF
C
C     Print and/or punch absolute cross sections
      WRITE (8,400)
400   FORMAT ('1Absolute cross sections to be fit')
      IF (IABSP.NE.0) WRITE (9,400)
      ADERCM=.FALSE.
      DO 907 I=1,NENERG
          CALL RDATA (I)
          ADERCM=ADERCM.OR.(ERCOM.NE.0.D0)
          IF (IABSP.NE.0) THEN
              WRITE (9,402) EX,NANG,NF,X0,OTEST,GNRM,ERCOM,ID
402           FORMAT (1X,F10.4,',',2(I3,','),2(1PE13.6,','),0PF10.6,
     +         ',',F10.5,',',A10)
              DO 908 J=1,NANG
                  WRITE (9,412) THX(J),SX(J),SXE(J),DELTH(J),
     +             DELPH(J),XMS(J)
412               FORMAT (1X,F9.4,5(',',1PE12.5))
908           CONTINUE
          ENDIF
          CALL PDATA (I,EX,NF,X0,OTEST,GNORM,0.D0,ERCOM,ID,THX,SX,SXE)
907   CONTINUE
      IF (IABSP.NE.0) WRITE (8,430)
430   FORMAT (/'0The above cross sections were punched')
C
C     Add in common error and reprint
      IF (ADERCM) THEN
          WRITE (8,400)
          WRITE (8,440)
440       FORMAT ('0A common (reproducability) error was added in'/
     +     ' quadrature to obtain the final cross section errors')
          DO 909 I=1,NENERG
              CALL RDATA (I)
              DO 910 J=1,NANG
                  SXE(J)=DSQRT(SXE(J)**2+(SX(J)*ERCOM)**2)
910           CONTINUE
              CALL PDATA (I,EX,NF,X0,OTEST,GNORM,0.D0,ERCOM,ID,THX,
     +         SX,SXE)
              CALL WDATA (I)
909       CONTINUE
      ENDIF
      IF (NGENM.NE.0) THEN
          WRITE (8,211)
211       FORMAT (/'0Generalized radial moment data'/
     +     '0Transition',T21,'ALPHA',T31,'Exponent',T43,'Moment',T57,
     +     'Error'/1X,60('-'))
          DO 911 I=1,NGENM
              WRITE (8,221) LEVNAM(I),GMOMAL(I),GMOMEX(I),RMOME(I),
     +         DRMOME(I)
221           FORMAT (1X,A,4F12.5)
911       CONTINUE
      ENDIF
      IF (INOFIT.NE.0) STOP 'FITR completed'
C
C     Initialize for fit
      WRITE (6,222)
222   FORMAT ('$Coulomb functions file? ')
      READ (5,11) FILE
      COULFI=FILE.NE.' '
      IF (COULFI) THEN
          INQUIRE (FILE=FILE,EXIST=RDCOUL)
          IF (RDCOUL) THEN
#             ifdef linux
C                 Linux g77 does not know about the "READONLY" attribute.
                  OPEN (12,FILE=FILE,STATUS='OLD',
     +             FORM='UNFORMATTED')
#             else
                  OPEN (12,FILE=FILE,STATUS='OLD',
     +             FORM='UNFORMATTED',READONLY)
#               endif
          ELSE
#             ifdef unix
                  OPEN (12,FILE=FILE,STATUS='UNKNOWN',
     +             FORM='UNFORMATTED')
#             else
C                 We assume VMS if not Unix.
                  OPEN (12,FILE=FILE,STATUS='NEW',
     +             FORM='UNFORMATTED')
#               endif
          ENDIF
      ENDIF
      DO 912 I=1,NENERG
          CALL RDATA (I)
          IF (RDCOUL) THEN
              READ (12,ERR=480,END=480) (FRN(J),FIN(J),GRN(J),
     +         GIN(J),J=1,NF)
              GO TO 481
480           STOP 'Invalid Coulomb functions input file'
          ELSE
              CALL COULS (Z0*ALPHA,X0,NF)
              IF (COULFI) WRITE (12) (FRN(J),FIN(J),GRN(J),GIN(J),J=1,NF)
          ENDIF
481       CALL WDATA (I)
912   CONTINUE
      RMS0=0.D0
      NITNRM=0
      MODELN=MODEL
      DO 913 I=1,NRAD
          RADIUS(I)=(I-1)*.1D0
913   CONTINUE
      DO 914 I=1,NPAR
          PARAM(I)=PAR(I)
914   CONTINUE
550   LAST=.FALSE.
      XLAM=.001D0
560   NIT=0
      CALL CPUBGN
C
C     Begin iterative fitting loop
570   NEWIT=.TRUE.
      IF (NIT.GE.NITMAX) LAST=.TRUE.
      LNGPRI=ILPR.NE.0.OR.(LAST.AND.INRM.EQ.0)
      IF (LNGPRI.OR.(NIT.EQ.0.AND.NITNRM.EQ.0)) THEN
C         Skip to new page on output file at start of fit and
C         after every iteration if a long print is being done.
          WRITE (8,11) '1'
      ELSE
          WRITE (8,11) ' '
      ENDIF
      WRITE (8,580) NIT
580   FORMAT (' Charge density iteration number',I3)
      NIT=NIT+1
      DO 915 I=1,NPAR
          PAR(I)=PARAM(I)
915   CONTINUE
      IF (LAST) XLAM=0.0D0
      DERIVB=DERIVB.OR.LAST
      WRITE (8,600) XLAM
600   FORMAT (' Proportion of gradient search included =',1PE7.0)
      IF (LNGPRI) WRITE (8,602)
602   FORMAT ('0',T4,'Group',T14,'Energy',T22,'THETA',T33,'q',T39,
     + 'q-eff',T47,'SIGMA',T58,'SIGMA',T74,'SIGMA',T90,'Difference',
     + T104,'Percent',T114,'CHI',T121,'CHI-',T128,'CHI-'/T56,
     + 'normalized',T72,'experimental',T105,'diff.',T119,'square',
     + T127,'square'/1X,131('-'))
      CHI=0.0D0
      NGR=0
      NGRPL=0
      QEXM=0.D0
      DO 916 I=1,NPARP2
          BETA(I)=0.D0
          DO 917 J=1,I
              CURV(I,J)=0.D0
917       CONTINUE
916   CONTINUE
C
C     Calculate curvature matrix (ALPHA) and vector (BETA)
C     First set up Lagrange multiplier constraint on Z.  Note that
C     the calculation of the constraint row and column is done
C     even if there is no constraint on the last iteration in
C     order to find the error in Z.
      IF (IZ.NE.0.OR.(IZ.EQ.0.AND.LAST.AND.NPARM1.GT.0)) THEN
          BETA(NPARP2)=Z0-Z
          DO 918 I=1,NPARM1
              IF (EPS(I).NE.0.0D0) THEN
                  DP=DPAR(I)*PAR(I)
                  IF (DP.EQ.0.D0) DP=.001D0
                  PARAM(I)=PAR(I)+DP
                  CURV(NPARP2,I)=(ZCALC(0)-Z)/DP
                  IF (DERIVB) THEN
C                     Calculate high-precision 2-point derivative
                      PARAM(I)=PAR(I)-DP
                      CURV(NPARP2,I)=((Z-ZCALC(0))/DP+
     +                 CURV(NPARP2,I))/2.D0
                  ENDIF
                  PARAM(I)=PAR(I)
              ENDIF
918       CONTINUE
      ENDIF
C
C     Calculate the portion of ALPHA and BETA for electron scattering
C     data and set up the elements for the RMS Lagrange multiplier
C     constraint.
      DO 919 I=1,NENERG
          CALL RDATA (I)
          IF (GNRM.NE.0) THEN
              IF (I.NE.1) THEN
C                 The following is done at the end of each normalization
C                 group except the last
                  WRITE (OUTPUT(124:131),690) CSQG
690               FORMAT (F8.2)
                  IF ((IDEVP.NE.0.OR.ICHP.NE.0).AND.LAST.AND.INRM.EQ.0)
     +             WRITE (9,700) IEND
700               FORMAT (T61,I10)
                  NGRPL=NGRPL+1
              ENDIF
              IF (GNRM.GT.0) THEN
                  NGR=NGR+1
                  SXT(NGR)=0.D0
                  STT(NGR)=0.D0
              ENDIF
              IF (LAST.AND.INRM.EQ.0) THEN
                  IF (IDEVP.NE.0) THEN
                      WRITE (9,720) TITLE,ID,TI,DA
720                   FORMAT ('   21~2Q~7EFF~8 (FM~5-1~6)',T67,'00E3',
     +                 T76,'000E3'/'   22% ~2DEVIATION FROM FIT',T67,
     +                 '00E5',T76,'000E5'/'FITR - ',A/
     +                 'P~2ERCENT DEVIATIONS FROM FIT FOR ~1',A/
     +                 T61,A,2X,A)
                      WRITE (9,700) IPOINT+NGRPL
                  ELSEIF (ICHP.NE.0) THEN
                      WRITE (9,721) TITLE,ID,TI,DA
721                   FORMAT ('   21~2Q~7EFF~8 (FM~5-1~6)',T67,'00E3',
     +                 T76,'000E3'/'    3~4C',T67,'00E4',T76,'000E4'/
     +                 'FITR - ',A/'D~2EVIATIONS FROM FIT IN UNITS ',
     +                 'OF ERROR BARS FOR ~1',A/T61,A,2X,A)
                      WRITE (9,700) IPOINT+NGRPL
                  ENDIF
              ENDIF
              CSQG=0.D0
          ENDIF
C         Write out the last line of the previous energy group.  This
C         line may contain the group CHI-square (see above). Then
C         initialize for the next energy group and perhaps normalization
C         group.
          IF (I.NE.1.AND.LNGPRI) THEN
              CALL NOSPAC(OUTPUT,J,K)
              WRITE (8,221) OUTPUT(1:K)
          ENDIF
          IF (GNRM.EQ.0) THEN
              OUTPUT=' '
          ELSE
              OUTPUT=ID
          ENDIF
          WRITE (OUTPUT(11:18),601) EX
601       FORMAT (F8.3)
          GAM=Z*ALPHA
          CALL RAVPRO (GAM,EX,EXCM,A,CONV,X0,STEP,OTEST,THX,RHO,
     +     ST,RHOERR,NANG,NF,NFL,DELTH,DELPH,XMS,0,RMS)
          CEFF=1.33333333333D0*GAM*HBARC/RMS
          IF (RMS0.EQ.0.D0) RMS0=RMS
          BETA(NPARP1)=RMS0-RMS
          IF (NPAR.GT.1) THEN
              DO 920 K=1,NPARM1
                  DO 921 J=1,NANG
                      DERIV(K,J)=0.0D0
921               CONTINUE
                  IF (EPS(K).NE.0.0D0) THEN
                      DP=DPAR(K)*PAR(K)
                      IF (DP.EQ.0.D0) DP=.001D0
                      PARAM(K)=PAR(K)+DP
                      CALL RAVPRO ((Z+CURV(NPARP2,K)*DP)*ALPHA,EX,
     +                 EXCM,A,CONV,X0,STEP,OTEST,THX,RHO,SR,RHOERR,
     +                 NANG,NF,NFL,DELTH,DELPH,XMS,0,RMSERR)
                      DO 922 J=1,NANG
                          DERIV(K,J)=GNORM*PAR(NPAR)*(SR(J)-ST(J))/DP
922                   CONTINUE
                      CURV(NPARP1,K)=(RMSERR-RMS)/DP
                      IF (DERIVB) THEN
                          PARAM(K)=PAR(K)-DP
                          CALL RAVPRO ((Z-CURV(NPARP2,K)*DP)*ALPHA,EX,
     +                     EXCM,A,CONV,X0,STEP,OTEST,THX,RHO,SR,RHOERR,
     +                     NANG,NF,NFL,DELTH,DELPH,XMS,0,RMSERR)
                          DO 923 J=1,NANG
                              DERIV(K,J)=(GNORM*PAR(NPAR)*
     +                         (ST(J)-SR(J))/DP+DERIV(K,J))/2.D0
923                       CONTINUE
                          CURV(NPARP1,K)=((RMS-RMSERR)/DP+
     +                     CURV(NPARP1,K))/2.D0
                      ENDIF
                      PARAM(K)=PAR(K)
                  ENDIF
920           CONTINUE
          ENDIF
          DO 924 J=1,NANG
              QQ=2.0D0*EX*DSIN(THX(J)*RAD/2.D0)/HBARC
              QEFF=QQ*(1.0D0+CEFF/EX)
              QEXM=DMAX1(QEXM,QEFF)
              SXCA=PAR(NPAR)*GNORM*ST(J)
              DIFY=SX(J)-SXCA
              CH=DIFY/SXE(J)
              CSQ=CH**2
              CHI=CHI+CSQ
              DIFYP=200.0D0*DIFY/(SX(J)+SXCA)
              SXEP=400.D0*SXE(J)*SXCA/(SX(J)+SXCA)**2
              IF (LAST) THEN
                  IF (GNORMS.GE.0.D0) THEN
                      STT(NGR)=STT(NGR)+(ST(J)*PAR(NPAR)/SXE(J))**2
                      SXT(NGR)=SXT(NGR)+SX(J)*ST(J)*PAR(NPAR)/SXE(J)**2
                  ENDIF
                  IF (INRM.EQ.0) THEN
                      IF (IDEVP.NE.0) THEN
                          WRITE (9,820) QEFF,DIFYP,SXEP
820                       FORMAT (E15.8,15X,2E15.8)
                      ELSEIF (ICHP.NE.0) THEN
                          WRITE (9,820) QEFF,CH,1.D0
                      ENDIF
                  ENDIF
              ENDIF
              IF (LNGPRI) THEN
                  CSQG=CSQG+CSQ
                  WRITE (OUTPUT(19:123),840) THX(J),QQ,QEFF,ST(J),
     +             SXCA,SX(J),SXE(J),DIFY,DIFYP,CH,CSQ
840               FORMAT (F8.3,2F8.4,1P3E11.4,'+-',1PE10.4,
     +             1PE12.4,0PF10.4,2F7.2)
                  IF (J.NE.NANG) THEN
C                     Write the line created above.  It may also contain
C                     the energy and group ID but not the group
C                     CHI-square.
                      CALL NOSPAC(OUTPUT,K,L)
                      WRITE (8,221) OUTPUT(1:L)
                      OUTPUT=' '
                  ENDIF
              ENDIF
              DERIV(NPAR,J)=GNORM*ST(J)
              DO 925 I1=1,NPAR
                  BETA(I1)=BETA(I1)+DERIV(I1,J)*(SX(J)-PAR(NPAR)*
     +             GNORM*ST(J))/SXE(J)**2
                  DO 926 J1=1,I1
                      CURV(I1,J1)=CURV(I1,J1)+DERIV(I1,J)*DERIV(J1,J)/
     +                 SXE(J)**2
926               CONTINUE
925           CONTINUE
924       CONTINUE
919   CONTINUE
C     Write out the final line for the last normalization group.
      IF (LNGPRI) THEN
          WRITE (OUTPUT(124:131),690) CSQG
          CALL NOSPAC(OUTPUT,I,J)
          WRITE (8,221) OUTPUT(1:J)
      ENDIF
      IF ((IDEVP.NE.0.OR.ICHP.NE.0).AND.LAST.AND.INRM.EQ.0)
     + WRITE (9,700) IEND
C
C     Calculate the portion of ALPHA and BETA for the generalized moment
C     data
      IF (NGENM.NE.0) THEN
          CSQG=0.
          IF (LNGPRI) WRITE (8,841)
841       FORMAT ('0Transition',T21,'ALPHA',T32,'Exponent',T46,'Moment',
     +     T75,'Moment',T90,'Difference',T104,'Percent',T114,'CHI',
     +     T121,'CHI-',T128,'CHI-'/T72,'experimental',T105,'diff.',
     +     T119,'square',T127,'square'/1X,131('-'))
          DO 927 I=1,NGENM
              RMOMT=GMOM(GMOMEX(I),GMOMAL(I))
              IF (NPARM1.GT.1) THEN
                  DO 928 J=1,NPARM1
                      DERIV(J,1)=0.0D0
                      IF (EPS(J).NE.0.0D0) THEN
                          DP=DPAR(J)*PAR(J)
                          IF (DP.EQ.0.D0) DP=.001D0
                          PARAM(J)=PAR(J)+DP
                          DERIV(J,1)=(GMOM(GMOMEX(I),GMOMAL(I))-
     +                     RMOMT)/DP
                          IF (DERIVB) THEN
                              PARAM(J)=PAR(J)-DP
                              DERIV(J,1)=((RMOMT-GMOM(GMOMEX(I),
     +                         GMOMAL(I)))/DP+DERIV(J,1))*.5D0
                          ENDIF
                          PARAM(J)=PAR(J)
                      ENDIF
928               CONTINUE
              ENDIF
              DIFY=RMOME(I)-RMOMT
              CH=DIFY/DRMOME(I)
              CSQ=CH**2
              CHI=CHI+CSQ
              IF (LNGPRI) THEN
                  OUTPUT=' '
                  DIFYP=200.0D0*DIFY/(RMOME(I)+RMOMT)
                  SXEP=400.D0*DRMOME(I)*RMOMT/(RMOME(I)+RMOMT)**2
                  CSQG=CSQG+CSQ
                  WRITE (OUTPUT(1:123),842) LEVNAM(I),GMOMAL(I),
     +             GMOMEX(I),RMOMT,RMOME(I),DRMOME(I),DIFY,DIFYP,CH,CSQ
842               FORMAT (A,F12.5,1X,F12.5,5X,1PE11.4,11X,1PE11.4,'+-',
     +             1PE10.4,1PE12.4,0PF10.4,2F7.2)
                  IF (I.EQ.NGENM) WRITE (OUTPUT(124:131),690) CSQG
                  CALL NOSPAC(OUTPUT,J,K)
                  WRITE (8,221) OUTPUT(1:K)
              ENDIF
              IF (NPARM1.GT.1) THEN
                  DO 929 I1=1,NPARM1
                      BETA(I1)=BETA(I1)+DERIV(I1,1)*(RMOME(I)-RMOMT)/
     +                 DRMOME(I)**2
                      DO 930 J1=1,I1
                          CURV(I1,J1)=CURV(I1,J1)+DERIV(I1,1)*
     +                     DERIV(J1,1)/DRMOME(I)**2
930                   CONTINUE
929               CONTINUE
              ENDIF
927       CONTINUE
      ENDIF
C
C     Calculate error matrix
860   I1=0
      DO 931 I=1,NPARP2
          IF (EPS(I).NE.0.0D0) THEN
              IF (CURV(I,I).LE.0.D0.AND.I.LE.NPAR) THEN
                  EPS(I)=0.0D0
                  WRITE (8,870) I
870               FORMAT ('0Charge density is insensitive to parameter '
     +             ,I2,'. It can''t be varied meaningfully.')
              ELSE
                  I1=I1+1
                  J1=0
                  DO 932 J=1,I
                      IF (EPS(J).NE.0.0D0) THEN
                          J1=J1+1
                          ALFA(I1,J1)=CURV(I,J)
                          ALFA(J1,I1)=ALFA(I1,J1)
                      ENDIF
932               CONTINUE
                  IF (I1.LE.NPAR) ALFA(I1,I1)=(1.0D0+XLAM)*ALFA(I1,I1)
              ENDIF
          ENDIF
931   CONTINUE
      IF (I1.GT.0) CALL INVRT2 (ALFA,I1,NPMXP3)
      I1=0
      DO 933 I=1,NPARP2
          DO 934 J=1,NPARP2
              ERROR(I,J)=0.0D0
934       CONTINUE
          IF (EPS(I).NE.0.0D0) THEN
              I1=I1+1
              J1=0
              DO 935 J=1,NPARP2
                  IF (EPS(J).NE.0.0D0) THEN
                      J1=J1+1
                      ERROR(I,J)=ALFA(I1,J1)
                  ENDIF
935           CONTINUE
          ENDIF
933   CONTINUE
      IF (NEWIT) THEN
C         Print out parameter information
          NEWIT=.FALSE.
          IF (LNGPRI) WRITE (8,11) ' '
          WRITE (8,940)
940       FORMAT (' Parameter',T21,'Value'/1X,60('-'))
          DO 936 I=1,NPAR
C             We use output to a character string here to avoid multiple
C             lines in the output file which end up being printed on
C             one line.
              OUTPUT=' '
              WRITE (OUTPUT,170) I,PAR(I)
              IF (EPS(I).NE.0.D0) THEN
                  OUTPUT(51:56)='varied'
                  IF (LAST) THEN
                      DEPS=DSQRT(ERROR(I,I))
                      WRITE (OUTPUT(30:47),970) DEPS
970                   FORMAT ('+-',1PE14.8)
                  ENDIF
              ELSE
                  OUTPUT(51:60)='NOT varied'
              ENDIF
              CALL NOSPAC(OUTPUT,J,K)
              WRITE (8,221) OUTPUT(1:K)
936       CONTINUE
          IF (LNGPRI) THEN
              WRITE (8,1010)
1010          FORMAT ('0Error matrix')
              DO 937 I=1,NPARP2,10
                  J=MIN0(I+9,NPARP2)
                  WRITE (8,1020) (K,K=I,J)
1020              FORMAT (4X,10(10X,I2))
                  WRITE (8,1025) ('------------',K=I,J)
1025              FORMAT (' -------',10A)
                  DO 938 K=1,NPARP2
                      WRITE (8,1035) K,(ERROR(K,L),L=I,J)
1035                  FORMAT (1X,I3,4X,10(1X,1PE11.4))
938               CONTINUE
                  WRITE (8,11) ' '
937           CONTINUE
              CALL CHGDE1 (1.D0,F1,F2,F3)
              DO 939 I=1,NRAD
                  CALL CHGDE2 (1.D0,RADIUS(I),RH)
                  RHO(I)=RH/(4.D0*PI)
939           CONTINUE
          ENDIF
          CALL CPUEND(TIMCPU,TIMCLK)
          CALL CPUBGN
          IF (LAST) GO TO 1200
          IF (LNGPRI) THEN
              WRITE (8,1050)
1050          FORMAT (1X,5('Radius',4X,'Charge',10X)/1X,123('-'))
              WRITE (8,1060) (RADIUS(I),RHO(I),I=1,NRAD)
1060          FORMAT (5(0PF6.1,1X,1PE13.6,6X))
          ENDIF
          IF (LNGPRI) WRITE (8,11) ' '
          WRITE (8,1070) CHI,Z,RMS,TIMCPU,TIMCLK
1070      FORMAT (' CHI-square = ',F10.3/' Z = ',F12.8,/
     +     ' RMS RADIUS = ',1PE13.7/' CPU time = ',0PF7.2,
     +     ' seconds'/' Clock time =',F7.2,' seconds')
      ENDIF
C
C     Step to new point in parameter space
      DO 969 I=1,NPAR
          PARAM(I)=PAR(I)
          DO 941 J=1,NPARP2
              PARAM(I)=PARAM(I)+ERROR(I,J)*BETA(J)
941       CONTINUE
969   CONTINUE
      ZERR=Z
      Z=ZCALC(0)
      CSQ=0.D0
      DO 942 I=1,NENERG
          CALL RDATA (I)
          CALL RAVPRO (Z*ALPHA,EX,EXCM,A,CONV,X0,STEP,OTEST,THX,
     +     RHO,ST,RHOERR,NANG,NF,NFL,DELTH,DELPH,XMS,0,RMS)
          DO 943 J=1,NANG
              CSQ=CSQ+((SX(J)-PARAM(NPAR)*GNORM*ST(J))/SXE(J))**2
943       CONTINUE
942   CONTINUE
      IF (NGENM.NE.0) THEN
          DO 944 I=1,NGENM
              CSQ=CSQ+((RMOME(I)-GMOM(GMOMEX(I),GMOMAL(I)))/
     +         DRMOME(I))**2
944       CONTINUE
      ENDIF
      IF (CSQ.LE.CHI) THEN
          XLAM=XLAM*.1D0
          LAMDWN=.TRUE.
          LAST=.TRUE.
          DO 945 I=1,NPAR
              LAST=LAST.AND.(DABS(PARAM(I)-PAR(I)).LE.
     +         DABS(EPS(I)*PAR(I)))
945       CONTINUE
      ELSE
          Z=ZERR
          IF (LAMDWN.AND..NOT.DERIVB) THEN
              WRITE (8,1130)
1130          FORMAT ('0Proportion of gradient search increased',
     +         ' after once decreasing.  From now on, 2-point',
     +         ' derivitives will be calculated and the'/
     +         ' preceding iteration will be redone.')
              DERIVB=.TRUE.
          ELSE
              WRITE (8,1120) CSQ
1120          FORMAT ('0Proportion of gradient search increased by a',
     +         ' factor of 10.  New CHI-square was',F10.3)
              XLAM=XLAM*10.D0
              IF (XLAM.LE.1.D0) GO TO 860
              WRITE (8,1170)
1170          FORMAT ('0Fitting is diverging!!',
     +         ' (try changing derivitive step)'/
     +         ' Preceding iteration will be re-done.')
              LAST=.TRUE.
          ENDIF
          DO 946 I=1,NPAR
              PARAM(I)=PAR(I)
946       CONTINUE
          NIT=NIT-1
      ENDIF
      GO TO 570
C
C     Calculate errors on derived quantities
1200  RMSERR=0.D0
      ZERR=0.D0
      DO 947 I=1,NRAD
          RHOERR(I)=0.D0
          IF (NPAR.GT.1) THEN
              DO 948 J=1,NPARM1
                  DERIV(J,I)=0.D0
                  IF (EPS(J).NE.0.D0) THEN
                      DP=DPAR(J)*PAR(J)
                      IF (DP.EQ.0.D0) DP=.001D0
                      PARAM(J)=PAR(J)+DP
                      CALL CHGDE1 (1.D0,F1,F2,F3)
                      CALL CHGDE2 (1.D0,RADIUS(I),DERIV(J,I))
                      PARAM(J)=PAR(J)-DP
                      CALL CHGDE1 (1.D0,F1,F2,F3)
                      CALL CHGDE2 (1.D0,RADIUS(I),RH)
                      DERIV(J,I)=(DERIV(J,I)-RH)/(8.D0*PI*DP)
                      PARAM(J)=PAR(J)
                  ENDIF
948           CONTINUE
          ENDIF
947   CONTINUE
      IF (NPAR.GT.1) THEN
          DO 949 I=1,NPARM1
              DO 950 J=1,NPARM1
                  RMSERR=RMSERR+ERROR(I,J)*CURV(NPARP1,J)*CURV(NPARP1,I)
                  ZERR=ZERR+ERROR(I,J)*CURV(NPARP2,J)*CURV(NPARP2,I)
                  DO 951 K=1,NRAD
                      RHOERR(K)=RHOERR(K)+ERROR(I,J)*
     +                 DERIV(J,K)*DERIV(I,K)
951               CONTINUE
950           CONTINUE
949       CONTINUE
          RMSERR=DSQRT(DDIM(RMSERR,0.D0))
          ZERR=DSQRT(DDIM(ZERR,0.D0))
          DO 952 I=1,NRAD
              RHOERR(I)=DSQRT(DDIM(RHOERR(I),0.D0))
952       CONTINUE
      ENDIF
      IF (LNGPRI) THEN
          WRITE (8,1260)
1260      FORMAT (1X,3('Radius',11X,'Charge',17X)/1X,113('-'))
          WRITE (8,1270) (RADIUS(I),RHO(I),RHOERR(I),I=1,NRAD)
1270      FORMAT (3(0PF6.1,1X,1PE13.6,'+-',1PE12.6,6X))
          WRITE (8,11) ' '
      ENDIF
      WRITE (8,1280) CHI,Z,ZERR,RMS,RMSERR,TIMCPU,TIMCLK
1280  FORMAT (' CHI-square = ',F10.3/' Z = ',F12.8,'+-',F9.8/
     + ' RMS radius = ',1PE13.7,'+-',1PE13.7/' CPU time = ',0PF7.2,
     + ' seconds'/' Clock time = ',F7.2,' seconds')
C
C     New normalization
      IF (INRM.NE.0) THEN
          NITNRM=NITNRM+1
          INRM=INRM-1
          IF (LNGPRI) THEN
              WRITE (8,400)
              WRITE (8,1480) NITNRM
1480          FORMAT (' Normalization iteration',I4)
          ELSE
              WRITE (8,11) ' '
              WRITE (8,1480) NITNRM
              WRITE (8,1482)
1482          FORMAT (' Group',T14,'Normalization'/1X,25('-'))
          ENDIF
          NGR=0
          DO 953 I=1,NENERG
              CALL RDATA (I)
              IF (GNRM.GT.0.D0) THEN
                  NGR=NGR+1
                  GNRM=SXT(NGR)/STT(NGR)
                  GNRMER=DSQRT(1.D0/STT(NGR))
                  GNORMS=GNRM
                  GNORM=GNORMS
              ELSE
                  GNRMER=0.D0
              ENDIF
              IF (LNGPRI) THEN
                  CALL PDATA (I,EX,NF,X0,OTEST,GNORM,GNRMER,
     +             ERCOM,ID,THX,SX,SXE)
              ELSE
                  WRITE (8,1481) ID,GNORM,GNRMER
1481              FORMAT (1X,A,F7.5,'+-',F6.5)
              ENDIF
              CALL WDATA (I)
953       CONTINUE
          IF (NGENM.NE.0.AND.LNGPRI) THEN
              WRITE (8,211)
              DO 954 I=1,NGENM
                  WRITE (8,221) LEVNAM(I),GMOMAL(I),GMOMEX(I),RMOME(I),
     +             DRMOME(I)
954           CONTINUE
          ENDIF
          IF (INRM.EQ.0) GO TO 560
          GO TO 550
      ENDIF
C
C     End of fit. Save the results that were ordered.
      IF (ICHP.NE.0.AND.IDEVP.NE.0) THEN
C         If the percent deviations were not plotted then if CHI plots
C         were requested they have already been done during the last
C         iteration. If both percent deviation and CHI plots are
C         requested then we must now do the CHI plots.  This requires
C         another pass of the data.
          NGRPL=0
          DO 955 I=1,NENERG
              CALL RDATA (I)
              IF (GNRM.NE.0) THEN
                  IF (I.NE.1) THEN
                      WRITE (9,700) IEND
                      NGRPL=NGRPL+1
                  ENDIF
                  WRITE (9,721) TITLE,ID,TI,DA
                  WRITE (9,700) IPOINT+NGRPL
              ENDIF
              GAM=Z*ALPHA
              CALL RAVPRO (GAM,EX,EXCM,A,CONV,X0,STEP,OTEST,THX,THR,
     +         ST,SR,NANG,NF,NFL,DELTH,DELPH,XMS,0,RMS)
              CEFF=1.33333333333D0*GAM*HBARC/RMS
              DO 956 J=1,NANG
                  QEFF=2.D0*EX*DSIN(.5D0*THX(J)*RAD)/HBARC*
     +             (1.D0+CEFF/EX)
                  CH=(SX(J)-PAR(NPAR)*GNORM*ST(J))/SXE(J)
                  WRITE (9,820) QEFF,CH,1.D0
956           CONTINUE
955       CONTINUE
          WRITE (9,700) IEND
      ENDIF
      IF (IRHOP.NE.0) THEN
C         Punch charge information
          WRITE (9,1290) TITLE,TI,DA
1290      FORMAT ('    8~2R (FM)',T67,'00E0',T76,'000E0'/
     +     '   26~4R~2(R) (PROTONS/FM~53~6)',T67,'00E1',T76,'000E1'/
     +     'FITR - ',A/
     +     'C~2HARGE DENSITY (SOLID) AND ERROR (DOTTED)'/
     +     T61,A,2X,A)
          WRITE (9,700) INULL
          DO 957 I=1,NRAD
              WRITE (9,820) RADIUS(I),RHO(I)
              IF (I.EQ.1) WRITE (9,700) ILINE
957       CONTINUE
          WRITE (9,700) INULL
          DO 958 I=1,NRAD
              WRITE (9,820) RADIUS(I),RHO(I)+RHOERR(I)
              IF (I.EQ.1) WRITE (9,700) IDOT
958       CONTINUE
          WRITE (9,700) INULL
          DO 959 I=1,NRAD
              WRITE (9,820) RADIUS(I),RHO(I)-RHOERR(I)
              IF (I.EQ.1) WRITE (9,700) IDOT
959       CONTINUE
          WRITE (9,700) IEND
          WRITE (9,1330) TITLE,TI,DA
1330      FORMAT ('    8~2R (FM)',T67,'00E0',T76,'000E0'/
     +     '   27~4DR~2(R) (PROTONS/FM~53~6)',T67,'00E2',T76,'000E2'/
     +     'FITR - ',A/
     +     'E~2RROR IN CHARGE DENSITY'/
     +     T61,A,2X,A)
          WRITE (9,700) INULL
          DO 960 I=1,NRAD
              WRITE (9,820) RADIUS(I),RHOERR(I)
              IF (I.EQ.1) WRITE (9,700) ILINE
960       CONTINUE
          WRITE (9,700) IEND
      ENDIF
      IF (IFFP.NE.0) THEN
C         Punch form factor info
          WRITE (9,1360) TITLE,TI,DA
1360      FORMAT ('   21~2Q~7EFF~8 (FM~5-1~6)',T67,'00E3',T76,'000E3'/
     +     '  -19~0M~2F(Q)~0M~2~52~6',T67,'00E4',T76,'000E4'/
     +     'FITR - ',A/
     +     'F~2ORM FACTOR SQUARED'/
     +     T61,A,2X,A)
          NGRPL=0
          IF (IFFP.NE.1) THEN
              EXM=HBARC*.75D0*QEXM-CEFF
              X0M=STEP*INT(3.0D0*RSQR*EXM/(HBARC*STEP))
              NFM=1.D0+3.0D0*X0M
              NFM=MIN(NFM,NFMX-1)
              OTESTM=180.0D0
              IF (2.D0*EXM*RSQR.GT.HBARC) OTESTM=2.0D0*RTD*
     +         DASIN(HBARC/(2.D0*EXM*RSQR))
          ENDIF
          CALL COULS (GAM,X0M,NFM)
          DO 961 I=1,NFM
              FRNR(I)=FRN(I)
              GRNR(I)=GRN(I)
              FINR(I)=FIN(I)
              GINR(I)=GIN(I)
961       CONTINUE
          DO 962 I=1,NENERG
              CALL RDATA (I)
              IF (EXM.NE.EX) THEN
                  CALL RAVPRO (GAM,EX,EXCM,A,CONV,X0,STEP,OTEST,THX,
     +             RHO,ST,RHOERR,NANG,NF,NFL,DELTH,DELPH,XMS,0,RMS)
                  DO 963 J=1,NANG
                      THX(J)=RTD*DACOS(1.D0-((EX+CEFF)/(EXM+CEFF))**2
     +                 *(1.D0-DCOS(THX(J)*RAD)))
                      DELTH(J)=0.D0
                      DELPH(J)=0.D0
                      XMS(J)=0.D0
963               CONTINUE
                  DO 964 J=1,NFM
                      FRN(J)=FRNR(J)
                      GRN(J)=GRNR(J)
                      FIN(J)=FINR(J)
                      GIN(J)=GINR(J)
964               CONTINUE
                  CALL RAVPRO (GAM,EXM,EXCM,A,CONV,X0M,STEP,OTESTM,THX,
     +             RHO,SR,RHOERR,NANG,NFM,NFL,DELTH,DELPH,XMS,0,RMS)
              ENDIF
              IF (GNRM.NE.0.D0) THEN
                  WRITE (9,700) IPOINT+NGRPL
                  NGRPL=NGRPL+1
              ENDIF
              DO 965 J=1,NANG
                  XLAM=(2.D0*EXM*DSIN(THX(J)*RAD/2.D0)**2/(GAM*
     +             HBARC*DCOS(THX(J)*RAD/2.D0)))**2/(GNORM*PAR(NPAR))
                  IF (EXM.NE.EX) XLAM=XLAM*SR(J)/ST(J)
                  QEFF=2.0D0*DSIN(THX(J)*RAD/2.D0)/HBARC*(EXM+CEFF)
                  SX(J)=SX(J)*XLAM
                  SXE(J)=SXE(J)*XLAM
                  WRITE (9,820) QEFF,SX(J),SXE(J)
965           CONTINUE
962       CONTINUE
          DO 966 I=1,100
              THX(I)=I*1.79D0
              DELTH(I)=0.D0
              DELPH(I)=0.D0
              XMS(I)=0.D0
966       CONTINUE
          DO 967 I=1,NFM
              FRN(I)=FRNR(I)
              GRN(I)=GRNR(I)
              FIN(I)=FINR(I)
              GIN(I)=GINR(I)
967       CONTINUE
          CALL RAVPRO (GAM,EXM,EXCM,A,CONV,X0M,STEP,OTESTM,THX,
     +     RHO,ST,RHOERR,100,NFM,NFL,DELTH,DELPH,XMS,0,RMS)
          WRITE (9,700) INULL
          DO 968 I=1,100
              QEFF=2.0D0*DSIN(THX(I)*RAD/2.D0)/HBARC*(EXM+CEFF)
              ST(I)=ST(I)*(2.D0*EXM*DSIN(THX(I)*RAD/2.D0)**2/
     +         (GAM*HBARC*DCOS(THX(I)*RAD/2.D0)))**2
              WRITE (9,820) QEFF,ST(I)
              IF (I.EQ.1) WRITE (9,700) ILINE
968       CONTINUE
          WRITE (9,700) IEND
      ENDIF
1475  STOP 'FITR completed'
      END
      SUBROUTINE PDATA (INIT,E,N,X,O,GN,GNE,EC,IDE,TH,S,SE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Include parameter file to define RAVPRO array sizes.
      INCLUDE 'ravpro.inc'
C     The maximum number of points in a given data set
      PARAMETER (NPTSMX=200)
      CHARACTER IDE*10,ID*10,IDS*10,OUTPUT*132
      DOUBLE PRECISION THX(NPTSMX),SX(NPTSMX),SXE(NPTSMX),DELTH(NPTSMX)
      DOUBLE PRECISION DELPH(NPTSMX),XMS(NPTSMX),TH(NPTSMX),S(NPTSMX)
      DOUBLE PRECISION SE(NPTSMX),FRN(NFMX),FIN(NFMX),GRN(NFMX),GIN(NFMX)
      COMMON /COUL/ FRN,GRN,FIN,GIN
      COMMON /DATA/ EX,NANG,NF,X0,OTEST,GNRM,GNORMS,GNORM,ERCOM,
     + ERCOMS,THX,SX,SXE,DELTH,DELPH,XMS,ID
      DATA NDIN/10/
C
C     Print the data
      IF (INIT.EQ.1) WRITE (8,10)
10    FORMAT ('0',T4,'Group',T14,'Normalization',T30,'Sys.',T35,
     + 'Energy',T43,'NF',T48,'X0',T54,'OTEST',T61,'Theta',T72,
     + 'Cross section',T92,'Percent',T100,'DELTA',T106,'DELTA',T114,
     + 'M. S.'/T30,'err.',T93,'error',T100,'THETA',T107,'PHI',T114,
     + 'data'/1X,121('-'))
      DO 969 J=1,NANG
          OUTPUT=' '
          IF (J.EQ.1.AND.GNRM.NE.0) THEN
              OUTPUT(1:10)=ID
              WRITE (OUTPUT(12:18),50) GN
50            FORMAT (F7.5)
              IF (GNE.NE.0) WRITE (OUTPUT(19:26),54) GNE
54            FORMAT ('+-',F6.5)
              WRITE (OUTPUT(27:32),55) EC
55            FORMAT (2PF5.2,'%')
          ENDIF
          IF (J.EQ.1) WRITE (OUTPUT(33:57),40) E,N,X,O
40        FORMAT (F8.3,I3,2F7.2)
          WRITE (OUTPUT(58:88),20) TH(J),S(J),SE(J)
20        FORMAT (F8.3,1PE11.4,'+-',1PE10.4)
          IF (S(J).NE.0.D0) WRITE (OUTPUT(89:97),30) SE(J)/S(J)
30        FORMAT (2PF8.2,'%')
          IF (DELTH(J).NE.0.D0) WRITE (OUTPUT(98:103),60) DELTH(J)
60        FORMAT (F6.3)
          IF (DELPH(J).NE.0.D0) WRITE (OUTPUT(104:109),60) DELPH(J)
          IF (XMS(J).GT.0.D0) WRITE (OUTPUT(110:119),80) XMS(J)
80        FORMAT (F6.3,' deg')
          IF (XMS(J).LT.0.D0) WRITE (OUTPUT(110:121),90) DABS(XMS(J))
90        FORMAT (F6.3,' g/cm2')
          WRITE (8,105) OUTPUT
105       FORMAT (1X,A)
969   CONTINUE
      RETURN
C
C     Read data from file
      ENTRY RDATA(INIT)
      IF (INIT.EQ.1) THEN
          NDIN=NDOUT
          REWIND NDIN
      ENDIF
      READ (NDIN) EX,NANG,NF,X0,OTEST,GNRM,ERCOM,ID
      IF (GNRM.NE.0) THEN
          GNORMS=GNRM
          ERCOMS=ERCOM
          IDS=ID
      ENDIF
      GNORM=DABS(GNORMS)
      ERCOM=ERCOMS
      ID=IDS
      READ (NDIN) (FRN(J),FIN(J),GRN(J),GIN(J),J=1,NF)
      READ (NDIN) (THX(J),SX(J),SXE(J),DELTH(J),DELPH(J),
     + XMS(J),J=1,NANG)
      RETURN
C
C     Write data to file
      ENTRY WDATA(INIT)
      IF (INIT.EQ.1) THEN
          NDOUT=21-NDIN
          REWIND NDOUT
      ENDIF
      WRITE (NDOUT) EX,NANG,NF,X0,OTEST,GNRM,ERCOM,ID
      WRITE (NDOUT) (FRN(J),FIN(J),GRN(J),GIN(J),J=1,NF)
      WRITE (NDOUT) (THX(J),SX(J),SXE(J),DELTH(J),DELPH(J),
     + XMS(J),J=1,NANG)
      RETURN
      END
      FUNCTION ZCALC(IDUM)
C     This routine uses the DEQ routines to find the integral of
C     the charge density from 0 to a cutoff. The single argument is
C     just a dummy so this can be a function.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RHOINT
      DOUBLE PRECISION R(3),YR(1),XINT(1),DEQTMP(1)
      COMMON /RHOI/ R,YR,GMOMEX,GMOMAL
      CALL CHGDE1 (1.D0,F1,F2,F3)
      R(1)=0.D0
      R(2)=.1D0
      R(3)=15.D0
      XINT(1)=0.D0
      CALL INDEQ2 (0)
      CALL DEQ2 (I,R,XINT,YR,DEQTMP,1,RHOINT)
      ZCALC=XINT(1)
      RETURN
      END
      SUBROUTINE RHOINT
C     This routine is used by the routine ZCALC via DEQ to
C     calculate the integral of the charge distribution from 0
C     to infinity.  CHGDE1 should already have been called.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION R(3),YR(1)
      COMMON /RHOI/ R,YR,GMOMEX,GMOMAL
      CALL CHGDE2 (1.D0,R(1),YR(1))
      YR(1)=YR(1)*R(1)*R(1)
      RETURN
      END
      FUNCTION GMOM(EX,AL)
C     This routine uses the DEQ routines to find the value of the
C     generalized moment defined by Ford and Rinker (Phys Rev. C7
C     (1973), 1208). Note that two integrals must be done, the integral
C     which gives the generalized moment and the integral of the charge.
C     This last is necessary because Z may change when charge density
C     parameters are varied.  EX is the exponent of r and AL is
C     the constant ALPHA in the exponential in the definition of a
C     generalized moment.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL GMOINT
      DOUBLE PRECISION R(3),YR(1),XINT(1),DEQTMP(1)
      COMMON /RHOI/ R,YR,GMOMEX,GMOMAL
C     First find the value of Z
      Z=ZCALC(0)
C     Now do the integral in the numerator of the definition
      R(1)=0.D0
      R(2)=.1D0
      R(3)=15.D0
      XINT(1)=0.D0
      GMOMEX=EX
      GMOMAL=AL
      CALL INDEQ2 (0)
      CALL DEQ2 (I,R,XINT,YR,DEQTMP,1,GMOINT)
      GMOM=(XINT(1)/Z)**(1./EX)
      RETURN
      END
      SUBROUTINE GMOINT
C     This routine is used by the routine GMOM via DEQ to
C     calculate the integral in the numerator of the definition of
C     the generalized moment.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION R(3),YR(1)
      COMMON /RHOI/ R,YR,GMOMEX,GMOMAL
      CALL CHGDE2 (1.D0,R(1),YR(1))
      YR(1)=YR(1)*R(1)**(2.+GMOMEX)*DEXP(-R(1)*GMOMAL)
      RETURN
      END
