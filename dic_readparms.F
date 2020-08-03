C $Header: /u/gcmpack/MITgcm/pkg/dic/dic_readparms.F,v 1.10 2010/04/11 20:54:49 jmc Exp $
C $Name:  $

#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: DIC_READPARMS
C !INTERFACE: ==========================================================
      SUBROUTINE DIC_READPARMS( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R DIC_READPARMS
C     | o Initialise and read dic package parameters
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIC_VARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid    :: My Thread Id. number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIC

C     === Local variables ===
C     msgBuf    :: Informational/error message buffer
C     iUnit     :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit, iter

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C-- Abiotic dic parameters:
C permil   :: set carbon mol/m3 <---> mol/kg conversion factor
C             default permil = 1024.5 kg/m3
C Pa2Atm   :: Conversion factor for atmospheric pressure pLoad (when coupled
C             to atmospheric model) into Atm. Default assumes pLoad in Pascal
C             1 Atm = 1.01325e5 Pa = 1013.25 mb

      NAMELIST /ABIOTIC_PARMS/ permil, Pa2Atm

#ifdef DIC_BIOTIC

C-- Biotic dic parameters:
C   DOPfraction :: fraction of new production going to DOP
C   KDOPRemin   :: DOP remineralization rate (1/s) = 1/(6 month)
C   KRemin      :: remin power law coeff
C   zcrit       :: Minimum Depth (m) over which biological activity
C                  is computed --> determines nlev as the indice of the
C                  first layer deeper than -zcrit
C   O2crit      :: critical oxygen level (mol/m3)
C   R_OP, R_CP  :: stochiometric ratios
C   R_NP, R_FeP, R_pop2poc
C   zca         :: scale depth for CaCO3 remineralization (m)
CC Parameters for light/nutrient limited bioac
C   parfrac     :: fraction of Qsw that is PAR
C   k0          :: light attentuation coefficient (1/m)
C   lit0        :: half saturation light constant (W/m2)
C   KPO4        :: half saturation phosphate constant (mol/m3)
C   KFE         :: half saturation fe constant (mol/m3)
CC Iron chemisty values
C   alpfe       :: solubility of aeolian fe
C   sedfe_factor :: ratio of sediment iron to sinking organic matter
C   sedfe_offset    :: y-axis crossing for Fe_flux = sedfe_factor*pflux + sedfe_offset
C   freefemax   :: max solubility of free iron (mol/m3)
CC Control variables
C   KScav       :: iron scavenging rate (s-1)
C   ligand_stab :: ligand-free iron stability constant (m3/mol)
C   ligand_tot  :: total free ligand  (mol/m3)
CC JML Added some extras for particle dependent Fe scavenging rates
CC See Parekh etal 2005, GBC equation 11
C   KScav_surf :: surface scavenging rate (s-1)
C   poc_wsp    :: Constant particulate sinking rate (m/s)
C   KScav_poc  :: ????
C   scav_exp    :: Honeyman exponent
C   alpha       :: timescale for biological activity
C                  read in alphaUniform and filled in 2d array alpha
C   rain_ratio  :: inorganic/organic carbon rain ratio
C                  read in rainRatioUniform and filled in 2d array rain_ratio

      NAMELIST /BIOTIC_PARMS/
     & DOPfraction, KDOPRemin, KRemin, zcrit, zca,
     & O2crit, R_OP, R_CP, R_NP, R_FeP, R_SIP,
     & parfrac, k0, lit0, KPO4, KFE, 
     & KScav, KScav_yr, KScav_surf, KScav_poc_yr, KScav_poc,
     & R_pop2poc, R_dust2fe, poc_wsp,  dust_wsp,
     & scav_exp, scav_ratio, KScav_dust_yr, 
     & KScav_dust, KScav_background, KScav_background_yr,
     & alpfe, freefemax, ligand_stab, ligand_tot, 
     & sedfe_factor, sedfe_offset, solfe, R_FeHe3,
     & fe_vent_depth_min, fe_sed_depth_max,
     & selectExportRemin,
#ifdef ALLOW_VARIABLE_LIGANDS     
     & gamma_lig, lambda_over_gamma, lambda_ref,
     & lig_q10, lig_tref,
#endif
     & alphaUniform_yr, alphaUniform, rainRatioUniform 
#endif

      NAMELIST /DIC_FORCING/
     &          DIC_windFile, DIC_atmospFile, DIC_iceFile,
     &          DIC_ironFile, DIC_silicaFile, 
#if ((defined ALLOW_FE) && (defined HYDROTHERMAL_FE))    
     &          DIC_hydroventFile,
#endif
     &          DIC_forcingPeriod, DIC_forcingCycle,
     &          dic_int1, dic_int2, dic_int3, dic_int4, dic_pCO2

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      _BEGIN_MASTER(myThid)

       permil      = 1. _d 0 / 1024.5 _d 0
       Pa2Atm      = 1.01325 _d 5

#ifdef DIC_BIOTIC
       DOPfraction = 0.67 _d 0

C Set export params default to the classic Martin curve
       selectExportRemin = 2      
C Set KRemin to UNSET_RL, it will be caught later and defaults set
C       KRemin(1)   = 0.9 _d 0
       DO iter=1,4
           KRemin(iter) = UNSET_RL
       ENDDO

       KDOPRemin   = 1. _d 0/(6. _d 0*30. _d 0*86400. _d 0)
       KRemin      = 0.9 _d 0
       zcrit       = 500. _d 0
       O2crit      = 4. _d -3
       R_OP        =-170. _d 0
       R_CP        = 117. _d 0
       R_NP        = 16. _d 0
       R_FeP       = 4.68 _d -4
       R_SIP       = 15. _d 0
       zca         = 3500. _d 0
       parfrac     = 0.4 _d 0
       k0          = 0.02 _d 0
       lit0        = 30. _d 0
       KPO4        = 5. _d -4
       KFE         = 1.2 _d -7
       alpfe       = 0.01 _d 0
       freefemax   = 3. _d -7
       rainRatioUniform = 7. _d -2

C JML general point:
C      1/yr values are easier to put in the data.dic file. Set those
C      default values here. Also set the 1/s values to "unset" here,
C      and either read them in from data.dic (backwards compatible) or
C      calculate them using default/read-in 1/yr values later. 
C      Alpha productivity rate default is 2e-3 molP/m3/yr 
C      (as high as 6e-3 molP/m3/yr in Parekh etal 2005, GBC eqn3)       
       alphaUniform_yr  = 2. _d -3
       alphaUniform     = UNSET_RL
       
#ifdef ALLOW_FE
C JML Added var for particle/dust dependent Fe scavenging
C      KScav is the uniform value when scavenging is constant rate
       KScav_yr       = 0.19 _d 0
       KScav          = UNSET_RL
C      KScav can be modified, but use constant value here for the surface
       KScav_surf     = UNSET_RL
C      KScav_poc is the initial partical dependent scavenging rate     
C      in l^scav_exp.mg^-scav_exp.d-1 (See Parekh etal 2005, GBC eqn11, ==0.079*360)
       KScav_poc_yr   = 0. _d 0
       KScav_poc      = UNSET_RL
C      R_pop2poc converts POP in mol P.m-3.s-1 to POC in mg C.l-1
       R_pop2poc      = R_CP * 12. _d 0
C      poc_wsp is the sinking speed of POC used for for particle dependent scavenging 
       poc_wsp        = 2900. _d 0/(360. _d 0*86400. _d 0)
C      scav_exp is the Honeyman exponent for particle dependent scavenging  (==0.58)
       scav_exp       = 0.58 _d 0
C      scav_ratio converts Honeyman's thorium scavenging rates to iron (==0.2) 
       scav_ratio     = 0.2 _d 0   
C      Dust concentration-dependent scavenging requires iron input to be converted back
C        to dust concentration (mol Fe/m2/s to g dust/m2/s)...
       R_dust2fe      = 3.5 _d -2/58.845 _d 0
C      ...and a sinking speed (see Aumont et al 2015, GMD)  
       dust_wsp       = 2. _d 0
C      and an initial dust scavenging rate (m3/g/s) (==150*360/86400*360)
       KScav_dust_yr  = 0. _d 0
       KScav_dust     = UNSET_RL
C      KScav_background is an additional rate added on top of the other values 
       KScav_background_yr = 0. _d 0 
       KScav_background    = UNSET_RL

#ifdef HYDROTHERMAL_FE
C These values for hydrothermal input of iron from Tagliabue et al (2014)       
       R_FeHe3     = 4.5 _d 8
       solfe       = 2. _d -3
       fe_vent_depth_min = 750. _d 0
#endif

#if ((defined SEDIMENT_FE_ELROD) || (defined SEDIMENT_FE_DALE))    
C sedimentary iron source 
       fe_sed_depth_max = 6000. _d 0
C    For Elrod et al (2004) use:
       sedfe_factor = 6.8 _d -4       
       sedfe_offset = -0.5 _d -6 / 86400. _d 0
C    For Dale et al (2015) use:
C       sedfe_factor = 170. _d 0 / (86400. _d 0 * 1. _d 6)      
C       sedfe_offset      = 0
#endif
       
C Ligand parameters for the iron cycle       
       ligand_stab = 1. _d 8
       ligand_tot  = 1. _d -6
       
#ifdef ALLOW_VARIABLE_LIGANDS 
C JML Added vars for prognostic Ligands 
       ligand_tot        = 0. _d 0     
       gamma_lig         = 0. _d 0
       lambda_over_gamma = 0. _d 0
C lambda_ref is calculated in dic_init_fixed
       lambda_ref  = 0. _d 0
       lig_q10     = 2. _d 0
       lig_tref    = 3.489 _d 0
#endif /* ALLOW_VARIABLE_LIGANDS */
#endif /* ALLOW_FE */
#endif /* DIC_BIOTIC */

       DIC_windFile  = ' '
       DIC_atmospFile= ' '
       DIC_iceFile   = ' '
       DIC_ironFile  = ' '
       DIC_silicaFile= ' '
#if ((defined ALLOW_FE) && (defined HYDROTHERMAL_FE))
C JML include iron source due to hydrothermal input into bottom layer            
       DIC_hydroventFile=' '       
#endif      
       dic_int1    = 0
       dic_int2    = 0
       dic_int3    = 0
       dic_int4    = 0
       dic_pCO2    = 278. _d -6

c default periodic forcing to same as for physics
       DIC_forcingPeriod = externForcingPeriod
       DIC_forcingCycle  = externForcingCycle

      WRITE(msgBuf,'(A)') ' DIC_READPARMS: opening data.dic'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

      CALL OPEN_COPY_DATA_FILE( 'data.dic', 'DIC_READPARMS',
     O                          iUnit, myThid )

C--   Read parameters from open data file:

C-    Abiotic parameters
      READ(UNIT=iUnit,NML=ABIOTIC_PARMS)

#ifdef DIC_BIOTIC
C-    Biotic parameters
      READ(UNIT=iUnit,NML=BIOTIC_PARMS)
#endif

C-    forcing filenames and parameters
      READ(UNIT=iUnit,NML=DIC_FORCING)

      WRITE(msgBuf,'(A)')
     &   ' DIC_READPARMS: finished reading data.dic'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

C--   Close the open data file
      CLOSE(iUnit)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C-    derive other parameters:

#ifdef DIC_BIOTIC
       QSW_underice = .FALSE.
#ifdef USE_QSW_UNDERICE
       QSW_underice = .TRUE.
#elif (defined (USE_QSW))
C if using Qsw and seaice, then ice fraction is already
C taken into account
       IF ( useSEAICE ) QSW_underice = .TRUE.
       IF ( useThSIce ) QSW_underice = .TRUE.
#endif
       IF ( alphaUniform .EQ. UNSET_RL ) THEN
          alphaUniform   = alphaUniform_yr/(360. _d 0 * 86400. _d 0)
       ENDIF
      
#ifdef ALLOW_FE
       IF ( KScav .EQ. UNSET_RL ) THEN
          KScav          = KScav_yr/(360. _d 0*86400. _d 0)
       ENDIF
       IF ( KScav_surf .EQ. UNSET_RL ) THEN
          KScav_surf     = KScav
       ENDIF
       IF ( KScav_poc .EQ. UNSET_RL ) THEN
          KScav_poc      = KScav_poc_yr/(360. _d 0*86400. _d 0)
       ENDIF
       IF ( KScav_dust .EQ. UNSET_RL ) THEN
          KScav_dust     = KScav_dust_yr/(360. _d 0*86400. _d 0)
       ENDIF
       IF ( KScav_background .EQ. UNSET_RL ) THEN
          KScav_background= KScav_background_yr/(360. _d 0*86400. _d 0)
       ENDIF
#endif /* ALLOW_FE */
#endif /* DIC_BIOTIC */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Print out parameter values :

      iUnit = standardMessageUnit
      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A)') '// ==================================='
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A)') '// DIC package parameters :'
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A)') '// ==================================='
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)

C- namelist ABIOTIC_PARMS
       CALL WRITE_0D_RL( permil, INDEX_NONE,'permil =',
     &  ' /* Ref. density to convert mol/m3 to mol/kg */')
       CALL WRITE_0D_RL( Pa2Atm, INDEX_NONE,'Pa2Atm =',
     &  ' /* Atmosph. pressure conversion coeff (to Atm) */')

#ifdef DIC_BIOTIC
C- namelist BIOTIC_PARMS
       CALL WRITE_0D_RL( DOPfraction, INDEX_NONE,'DOPfraction =',
     &  ' /* Fraction of new production going to DOP */')
       CALL WRITE_0D_RL( KDOPRemin, INDEX_NONE,'KDOPRemin =',
     &  ' /* DOP remineralization rate (1/s) */')

       IF ( KRemin(1) .EQ. UNSET_RL ) THEN
C Set some default values depending on the model chosen (default 2)
          IF ( selectExportRemin .EQ. 1 ) THEN
C Exponential model (EQN1 Cael & Bisson) minimizing absolute error
             KRemin(1)=1.518 _d 0
             KRemin(2)=2.394 _d 0
             KRemin(3)=0. _d 0
             KRemin(4)=0. _d 0
          ELSEIF (selectExportRemin .EQ. 2) THEN
C Power Law model/ Martin curve (EQN2 Cael & Bisson) 
             KRemin(1)=0.9 _d 0
          ELSEIF ( selectExportRemin .EQ. 3 ) THEN
C Ballast model (EQN3 Cael & Bisson) minimizing absolute error
             KRemin(1)=1.559 _d 0
             KRemin(2)=2.006 _d 0
             KRemin(3)=0.03394 _d 0
             KRemin(4)=0. _d 0
          ELSEIF ( selectExportRemin .EQ. 4 ) THEN
C Rational model (EQN4 Cael & Bisson) minimizing absolute error
             KRemin(1)=1.2228 _d 0
             KRemin(2)=0.2228 _d 0
             KRemin(3)=0. _d 0
             KRemin(4)=0. _d 0
          ELSEIF ( selectExportRemin .EQ. 5 ) THEN
C Double exponential model (EQN5 Cael & Bisson) minimizing absolute error
             KRemin(1)=1.669 _d 0
             KRemin(2)=1.3372 _d 0
             KRemin(3)=0.1323 _d 0
             KRemin(4)=19.42 _d 0
          ELSEIF ( selectExportRemin .EQ. 6 ) THEN
C Stretched exponential (EQN6 Cael & Bisson) minimizing absolute error
             KRemin(1)=7387 _d 0
             KRemin(2)=8.934719866289981 _d 0
             KRemin(3)=1.0 _d 0
             KRemin(4)=0. _d 0
          ELSEIF ( selectExportRemin .EQ. 7 ) THEN
C Gamma model (EQN7 Cael & Bisson) minimizing absolute error
             KRemin(1)=0.6700 _d 0
             KRemin(2)=6.888 _d 0
             KRemin(3)=0. _d 0
             KRemin(4)=0. _d 0
          ENDIF
       ENDIF
       
       IF ( selectExportRemin .EQ. 1 ) THEN
C Exponential model (EQN1 Cael & Bisson)
           CALL WRITE_0D_RL( KRemin(1), INDEX_NONE,'KRemin(1) =',
     &  ' /* Remin exponential model coeff. */')
           CALL WRITE_0D_RL( KRemin(2), INDEX_NONE,'KRemin(2) =',
     &  ' /* Remin exponential model coeff. */')
       ELSEIF ( selectExportRemin .EQ. 2 ) THEN
C Power Law model/ Martin curve (EQN2 Cael & Bisson) 
           CALL WRITE_0D_RL( KRemin(1), INDEX_NONE,'KRemin =',
     &  ' /* Remin power law coeff. */')
       ELSEIF ( selectExportRemin .EQ. 3 ) THEN
C Ballast model (EQN3 Cael & Bisson)
           CALL WRITE_0D_RL( KRemin(1), INDEX_NONE,'KRemin(1) =',
     &  ' /* Remin ballast model coeff. */')
           CALL WRITE_0D_RL( KRemin(2), INDEX_NONE,'KRemin(2) =',
     &  ' /* Remin ballast model coeff. */')
           CALL WRITE_0D_RL( KRemin(3), INDEX_NONE,'KRemin(3) =',
     &  ' /* Remin ballast model coeff. */')
       ELSEIF ( selectExportRemin .EQ. 4 ) THEN
C Rational model (EQN4 Cael & Bisson)
           CALL WRITE_0D_RL( KRemin(1), INDEX_NONE,'KRemin(1) =',
     &  ' /* Remin rational model coeff. */')
           CALL WRITE_0D_RL( KRemin(2), INDEX_NONE,'KRemin(2) =',
     &  ' /* Remin rational model coeff. */')
       ELSEIF ( selectExportRemin .EQ. 5 ) THEN
C Double exponential model (EQN5 Cael & Bisson)
           CALL WRITE_0D_RL( KRemin(1), INDEX_NONE,'KRemin(1) =',
     &  ' /* Remin double exponential model coeff. */')
           CALL WRITE_0D_RL( KRemin(2), INDEX_NONE,'KRemin(2) =',
     &  ' /* Remin double exponential coeff. */')
           CALL WRITE_0D_RL( KRemin(3), INDEX_NONE,'KRemin(3) =',
     &  ' /* Remin double exponential coeff. */')
           CALL WRITE_0D_RL( KRemin(4), INDEX_NONE,'KRemin(4) =',
     &  ' /* Remin double exponential coeff. */')
       ELSEIF ( selectExportRemin .EQ. 6  ) THEN
C Stretched exponential (EQN6 Cael & Bisson)
           CALL WRITE_0D_RL( KRemin(1), INDEX_NONE,'KRemin(1) =',
     &  ' /* Remin stretched exponential model coeff. */')
           CALL WRITE_0D_RL( KRemin(2), INDEX_NONE,'KRemin(2) =',
     &  ' /* Remin stretched exponential coeff. */')
           CALL WRITE_0D_RL( KRemin(3), INDEX_NONE,'KRemin(3) =',
     &  ' /* Remin stretched exponential coeff. */')
       ELSEIF ( selectExportRemin .EQ. 7 ) THEN
C Gamma model (EQN7 Cael & Bisson)
           CALL WRITE_0D_RL( KRemin(1), INDEX_NONE,'KRemin(1) =',
     &  ' /* Remin gamma model coeff. */')
           CALL WRITE_0D_RL( KRemin(2), INDEX_NONE,'KRemin(2) =',
     &  ' /* Remin gamma model coeff. */')
       ELSE
C Invalid option selected, raise an error        
           WRITE(msgBuf,'(A)')
     &          'DIC_READPARMS: selectExportRemin 0< >7'
           CALL PRINT_ERROR( msgBuf, myThid )
           WRITE(msgBuf,'(A)')
     &          'DIC_READPARMS: valid options are 1-7'
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: selectExportRemin error' 
       ENDIF

C       CALL WRITE_0D_RL( KRemin, INDEX_NONE,'KRemin =',
C     &  ' /* Remin power law coeff. */')
       CALL WRITE_0D_RL( zcrit, INDEX_NONE,'zcrit =',
     &  ' /* Minimum depth for biological activity (m) */')
       CALL WRITE_0D_RL( O2crit, INDEX_NONE,'O2crit =',
     &  ' /* Critical oxygen level (mol/m3) */')
       CALL WRITE_0D_RL( R_OP, INDEX_NONE,'R_OP =',
     &  ' /* Stochiometric ratio R_OP */')
       CALL WRITE_0D_RL( R_CP, INDEX_NONE,'R_CP =',
     &  ' /* Stochiometric ratio R_CP */')
       CALL WRITE_0D_RL( R_NP, INDEX_NONE,'R_NP =',
     &  ' /* Stochiometric ratio R_NP */')
       CALL WRITE_0D_RL( R_FeP, INDEX_NONE,'R_FeP =',
     &  ' /* Stochiometric ratio R_FeP */')
       CALL WRITE_0D_RL( zca, INDEX_NONE,'zca =',
     &  ' /* Scale depth for CaCO3 remineralization (m) */')
       CALL WRITE_0D_RL( parfrac, INDEX_NONE,'parfrac =',
     &  ' /* Fraction of Qsw that is PAR */')
       CALL WRITE_0D_RL( k0, INDEX_NONE,'k0 =',
     &  ' /* Light attentuation coefficient (1/m) */')
       CALL WRITE_0D_RL( lit0, INDEX_NONE,'lit0 =',
     &  ' /* Half saturation light constant (W/m2) */')
       CALL WRITE_0D_RL( KPO4, INDEX_NONE,'KPO4 =',
     &  ' /* Half saturation phosphate constant (mol/m3) */')
#ifdef ALLOW_FE     
       CALL WRITE_0D_RL( KFE, INDEX_NONE,'KFE =',
     &  ' /* Half saturation fe constant (mol/m3) */')
       CALL WRITE_0D_RL( alpfe, INDEX_NONE,'alpfe =',
     &  ' /* Solubility of aeolian fe */')
       CALL WRITE_0D_RL( freefemax, INDEX_NONE,'freefemax =',
     &  ' /* Max solubility of free iron (mol/m3) */')

#ifdef ALLOW_VARIABLE_LIGANDS
C JML Prognostic Ligand parameters
      IF ( gamma_lig.NE.0. _d 0 ) THEN
            CALL WRITE_0D_RL(gamma_lig, INDEX_NONE,'gamma_lig =',
     &  ' /* Fraction of biomass that produces ligands */')   
            CALL WRITE_0D_RL(lambda_over_gamma, INDEX_NONE,
     &       'lambda_over_gamma =',
     &  ' /* Lambda over Gamma, ratio of ligand coefficients */')      
             CALL WRITE_0D_RL(lambda_over_gamma*gamma_lig, 
     &        INDEX_NONE,'lambda_ref =',
     &  ' /* reference ligand degradation timescale */')   
            CALL WRITE_0D_RL(lig_q10, INDEX_NONE,'lig_q10 =',
     &  ' /* Variable ligand degradation Q10 value */')   
            CALL WRITE_0D_RL(lig_tref, INDEX_NONE,'lig_tref =',
     &  ' /* Variable ligand degradation reference temperature */')   
      ENDIF
#endif /* ALLOW_VARIABLE_LIGANDS */

#ifdef HYDROTHERMAL_FE
C JML Hydrothermal iron input parameters
        IF ( DIC_hydroventFile .NE. ' '  ) THEN
            CALL WRITE_0D_RL(R_FeHe3, INDEX_NONE,'R_FeHe3 =',
     &  ' /* Ratio of iron to 3He in vent fluid */')   
            CALL WRITE_0D_RL(solfe, INDEX_NONE,'solfe =',
     &  ' /* Solubility of hydrothermal iron */')   
            CALL WRITE_0D_RL(fe_vent_depth_min, 
     &                INDEX_NONE,'fe_vent_depth_min =',
     &  ' /* Minimum depth that hydrothermal iron is input */')   
      ENDIF
#endif /* HYDROTHERMAL_FE */

#if ((defined SEDIMENT_FE_ELROD) || (defined SEDIMENT_FE_DALE))
C JML Sediment iron input parameters
            CALL WRITE_0D_RL(sedfe_factor, INDEX_NONE,'sedfe_factor =',
     &  ' /* Gradient of sediment iron-poc relationship */')   
            CALL WRITE_0D_RL(sedfe_offset, INDEX_NONE,'sedfe_offset =',
     &  ' /* Intercept of sediment iron-poc relationship */')   
            CALL WRITE_0D_RL(fe_sed_depth_max, INDEX_NONE,
     &  'fe_sed_depth_max =',
     &  ' /*  Maximum depth that sediment iron is input */')   
#endif /* SEDIMENT_FE */

C JML Fe Scavenging parameters
       CALL WRITE_0D_RL( KScav_yr, INDEX_NONE,'KScav_yr =',
     &  ' /* Iron scavenging rate (1/yr)*/')
       CALL WRITE_0D_RL( KScav_background_yr, INDEX_NONE,
     &  'KScav_background_yr =',
     &  ' /* Background scavenging rate (1/yr)*/')
       CALL WRITE_0D_RL( KScav_dust_yr, INDEX_NONE,'KScav_dust_yr =',
     &  ' /* Iron scavenging rate due to dust (1/yr)*/')
       CALL WRITE_0D_RL( dust_wsp, INDEX_NONE,'dust_wsp',
     &  ' /* Dust particle sinking speed (m/s) */')
       CALL WRITE_0D_RL( R_dust2fe, INDEX_NONE,'R_dust2fe =',
     &  ' /* Content of iron in dust (3.5%/58gFe/mol) */')
C JML write out particle depended scavenging parameters
      IF ( KScav_poc_yr.GT.0. _d 0 ) THEN
            CALL WRITE_0D_RL(KScav_poc_yr, INDEX_NONE,'KScav_poc_yr =',
     &  ' /* Iron scavenging rate due to organic particles (1/yr)*/')
            CALL WRITE_0D_RL( scav_exp, INDEX_NONE,'scav_exp =',
     &  ' /* Honeyman exponent for scavenging */')
            CALL WRITE_0D_RL( scav_ratio, INDEX_NONE,'scav_ratio =',
     &  ' /* Honeyman ratio of thorium to iron scavenging */')
            CALL WRITE_0D_RL( poc_wsp, INDEX_NONE,'poc_wsp',
     &  ' /* Organic particle sinking speed (m/s) */')
            CALL WRITE_0D_RL( R_pop2poc, INDEX_NONE,'R_pop2poc =',
     &  ' /* Stochiometric ratio R_pop2poc */')
      ENDIF
       CALL WRITE_0D_RL( ligand_stab, INDEX_NONE,'ligand_stab =',
     &  ' /* Ligand-free iron stability constant (m3/mol) */')
       CALL WRITE_0D_RL( ligand_tot, INDEX_NONE,'ligand_tot =',
     &  ' /* Total free ligand  (mol/m3) */')
#endif /* ALLOW_FE */       
       CALL WRITE_0D_RL( alphaUniform, INDEX_NONE,'alphaUniform =',
     &  ' /* Timescale for biological activity */')
       CALL WRITE_0D_RL(rainRatioUniform,INDEX_NONE,'rainRatioUniform=',
     &  ' /* Inorganic/organic carbon rain ratio */')

       CALL WRITE_0D_L( QSW_underice, INDEX_NONE, 'QSW_underice  =',
     &  '  /* Flag for Qsw under Sea-Ice (i.e. SI fract included) */')
#endif

C- namelist DIC_FORCING
       CALL WRITE_0D_C( DIC_windFile, -1, INDEX_NONE, 'DIC_windFile =',
     & '  /* File name of wind speeds */')
       CALL WRITE_0D_C( DIC_atmospFile, -1,INDEX_NONE,'DIC_atmospFile=',
     & '  /* File name of atmospheric pressure*/')
       CALL WRITE_0D_C( DIC_iceFile, -1, INDEX_NONE, 'DIC_iceFile =',
     & '  /* File name of seaice fraction */')
#ifdef ALLOW_FE
       CALL WRITE_0D_C( DIC_ironFile, -1, INDEX_NONE, 'DIC_ironFile =',
     & '  /* File name of aeolian iron flux */')
#ifdef HYDROTHERMAL_FE
       CALL WRITE_0D_C( DIC_hydroventFile, -1,INDEX_NONE,
     & 'DIC_hydroventFile=',
     & '  /* File name of Hydrothermal 3He input */')
#endif /* HYDROTHERMAL_FE */
#endif /* ALLOW_FE */
       CALL WRITE_0D_C( DIC_silicaFile, -1,INDEX_NONE,'DIC_silicaFile=',
     & '  /* File name of surface silica */')
       CALL WRITE_0D_RL( DIC_forcingPeriod,
     &   INDEX_NONE,'DIC_forcingPeriod =',
     &  ' /* Periodic forcing parameter specific for DIC (s) */')
       CALL WRITE_0D_RL( DIC_forcingCycle,
     &   INDEX_NONE,'DIC_forcingCycle =',
     &  ' /* Periodic forcing parameter specific for DIC (s) */')
       CALL WRITE_0D_I( dic_int1, INDEX_NONE, 'dic_int1 =',
     &  '  /*  */')
       CALL WRITE_0D_I( dic_int2, INDEX_NONE, 'dic_int2 =',
     &  '  /*  */')
       CALL WRITE_0D_I( dic_int3, INDEX_NONE, 'dic_int3 =',
     &  '  /*  */')
       CALL WRITE_0D_I( dic_int4, INDEX_NONE, 'dic_int4 =',
     &  '  /*  */')
       CALL WRITE_0D_RL( dic_pCO2, INDEX_NONE,'dic_pCO2 =',
     &  ' /* Atmospheric pCO2 to be read in data.dic */')

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( dic_int1.EQ.0 .AND. dic_pCO2.NE.278. _d -6 ) THEN
        WRITE(msgBuf,'(A)')
     &    'DIC_READPARMS: cannot change default dic_pCO2 if dic_int1=0'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R DIC_READPARMS: dic_pCO2 error'
      ENDIF

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#endif /* ALLOW_DIC */

      RETURN
      END
