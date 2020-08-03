C $Header: /u/gcmpack/MITgcm/pkg/dic/DIC_VARS.h,v 1.8 2011/04/19 21:34:32 stephd Exp $
C $Name:  $

#include "DIC_OPTIONS.h"

C     *==========================================================*
C     | DIC_VARS.h
C     | o Carbon Variables
C     *==========================================================*

       COMMON /CARBON_NEEDS/
     &              AtmospCO2, AtmosP, pH, pCO2,
     &              co3, hco3, co2aq, FluxCO2,
     &              wind, FIce, Silica, Kwexch_Pre
      _RL  AtmospCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  AtmosP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  pH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  pCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  co3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  hco3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  co2aq(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FluxCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  wind(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Kwexch_Pre(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Silica(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      
       COMMON /CARBON_CHEM/
     &                     ak0,ak1,ak2,akw,akb,aks,akf,
     &                     ak1p,ak2p,ak3p,aksi, fugf, 
     &                     ff,ft,st,bt, Ksp_TP_Calc
      _RL  ak0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  akw(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  akb(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  aks(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  akf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak1p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak2p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak3p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  aksi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C Fugacity Factor added by Val Bennington Nov. 2010
      _RL  fugf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ft(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  st(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  bt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Ksp_TP_Calc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

       COMMON /OXYGEN_CHEM/
     &              oA0,oA1,oA2,oA3,oA4,oA5,
     &              oB0,oB1,oB2,oB3,
     &              oC0
      _RL oA0,oA1,oA2,oA3,oA4,oA5
      _RL oB0,oB1,oB2,oB3
      _RL oC0

C permil : is conversion factor for mol/m3 to mol/kg
C          assumes uniform (surface) density
C Pa2Atm : for conversion of atmospheric pressure
C          when coming from atmospheric model
       COMMON /GLOBAL_SURF_MEAN/
     &                          gsm_alk,gsm_s,gsm_t,gsm_dic,
     &                          gsm_c14,permil,Pa2Atm
      _RL  gsm_alk
      _RL  gsm_s
      _RL  gsm_t
      _RL  gsm_DIC
      _RL  gsm_C14
      _RL  permil
      _RL  Pa2Atm

C schmidt number coefficients
      COMMON /DIC_SCHMIDT_NO/
     &                    sca1, sca2, sca3, sca4,
     &                    sox1, sox2, sox3, sox4
      _RL  sca1
      _RL  sca2
      _RL  sca3
      _RL  sca4
      _RL  sox1
      _RL  sox2
      _RL  sox3
      _RL  sox4

C--   COMMON /DIC_FILENAMES/
C  DIC_windFile    :: file name of wind speeds
C  DIC_atmospFile  :: file name of atmospheric pressure
C  DIC_iceFile     :: file name of seaice fraction
C  DIC_ironFile    :: file name of aeolian iron flux
C  DIC_silicaFile  :: file name of surface silica
C  DIC_forcingPeriod :: periodic forcing parameter specific for dic (seconds)
C  DIC_forcingCycle  :: periodic forcing parameter specific for dic (seconds)
C  dic_pCO2          :: Atmospheric pCO2 to be rad in data.dic
C  dic_int*          :: place holder to read in a integer number, set at run time
C JML include iron source due to hydrothermal input into bottom layer            

      COMMON /DIC_FILENAMES/
     &        DIC_windFile, DIC_atmospFile, DIC_iceFile,
     &        DIC_ironFile, DIC_hydroventFile, DIC_silicaFile,
     &        DIC_forcingPeriod, DIC_forcingCycle,
     &        dic_int1, dic_int2, dic_int3, dic_int4, 
     &        dic_pCO2

      CHARACTER*(MAX_LEN_FNAM) DIC_windFile
      CHARACTER*(MAX_LEN_FNAM) DIC_atmospFile
      CHARACTER*(MAX_LEN_FNAM) DIC_iceFile
      CHARACTER*(MAX_LEN_FNAM) DIC_ironFile
      CHARACTER*(MAX_LEN_FNAM) DIC_hydroventFile
      CHARACTER*(MAX_LEN_FNAM) DIC_silicaFile
      _RL     DIC_forcingPeriod
      _RL     DIC_forcingCycle
      _RL dic_pCO2
      INTEGER dic_int1
      INTEGER dic_int2
      INTEGER dic_int3
      INTEGER dic_int4

#ifdef DIC_BIOTIC
C     *==========================================================*
C     | o Biological Carbon Variables
C     *==========================================================*
C JML Added some extras for particle dependent Fe scavenging rates
C       and prognostic Ligands
      COMMON /BIOTIC_NEEDS/
     &     BIOave, CARave, SURave, SUROave, pCO2ave, pHave,
     &     fluxCO2ave, omegaCave, pfluxave, epfluxave, cfluxave,
     &     DIC_timeAve,
     &     alpha, rain_ratio, InputFe, omegaC,
     &     Kpo4, DOPfraction, zcrit, KRemin, O2crit,
     &     KDOPremin,zca,R_op,R_cp,R_NP, R_FeP, R_SIP,
     &     alpfe, ligand_stab, ligand_tot, KFE, freefemax, 
     &     KScav, KScav_yr, KScav_surf, KScav_poc_yr, KScav_poc,
     &     R_pop2poc, R_dust2fe, poc_wsp,  dust_wsp,
     &     scav_exp, scav_ratio, KScav_dust_yr, 
     &     KScav_dust, KScav_background, KScav_background_yr,
     &     sedfe_factor, sedfe_offset, fe_sed_depth_max,
     &     HydroInputHe3, solfe, R_FeHe3, fe_vent_depth_min,
     &     par, parfrac, k0, lit0,
     &     alphaUniform_yr, alphaUniform, rainRatioUniform,
     &     alphamax, alphamin, 
     &     calpha, crain_ratio, cInputFe, calpfe, feload, cfeload,
     &     nlev, QSW_underice,
     &     selectExportRemin

      INTEGER selectExportRemin
      INTEGER nlev

C     For averages
      _RL BIOave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL CARave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL SURave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL SUROave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pCO2ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pHave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL fluxCO2ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL OmegaCave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL pfluxave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL epfluxave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL cfluxave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL DIC_timeAve(nSx,nSy)

C     values for biogeochemistry
C JML Added some extras for particle dependent Fe scavenging rates,
C       prognostic Ligands, and hajoon's sediment source
C   sedfe_factor :: ratio of sediment iron to sinking organic matter
C   sedfe_offset      :: Sediment Fe flux, intersect value in:
C                    Fe_flux = sedfe_factor*pflux + sedfe_offset
      _RL par(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL alpha(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL rain_ratio(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL InputFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL HydroInputHe3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL omegaC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL Kpo4
      _RL DOPfraction
      _RL zcrit
      _RL KRemin(4)
      _RL KDOPremin
      _RL zca
      _RL R_op
      _RL R_cp
      _RL R_NP
      _RL R_FeP
      _RL R_SIP
      _RL O2crit
      _RL alpfe
      _RL ligand_stab
      _RL ligand_tot
      _RL KFe
      _RL freefemax
C Scavenging values
      _RL KScav_yr
      _RL KScav
      _RL KScav_surf
      _RL KScav_poc_yr
      _RL KScav_poc
      _RL R_pop2poc
      _RL R_dust2fe
      _RL poc_wsp
      _RL dust_wsp
      _RL scav_exp
      _RL scav_ratio
      _RL KScav_dust_yr
      _RL KScav_dust
      _RL KScav_background_yr
      _RL KScav_background
      _RL sedfe_factor
      _RL sedfe_offset
      _RL R_FeHe3
      _RL solfe
      _RL fe_vent_depth_min
      _RL fe_sed_depth_max
C     values for light limited bio activity
      _RL k0, parfrac, lit0
      _RL alphaUniform_yr
      _RL alphaUniform
      _RL rainRatioUniform
      _RL alphamax, alphamin
      _RL calpha
      _RL crain_ratio
      _RL cInputFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL calpfe
      _RL cfeload
      _RL feload(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      LOGICAL QSW_underice

#ifdef ALLOW_VARIABLE_LIGANDS
      COMMON /PROGNOSTIC_LIGANDS/
     &     gamma_lig, lambda_over_gamma, lambda_ref,
     &     lig_q10, lig_tref
     
      _RL gamma_lig
      _RL lambda_over_gamma
      _RL lambda_ref
      _RL lig_q10
      _RL lig_tref
#endif
      
#endif /* DIC_BIOTIC */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
