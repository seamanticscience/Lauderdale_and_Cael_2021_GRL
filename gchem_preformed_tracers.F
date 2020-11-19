C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_forcing_sep.F,v 1.30 2011/12/22 19:05:27 jmc Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#ifdef ALLOW_DIC
#include "DIC_OPTIONS.h"
#endif
#ifdef ALLOW_DARWIN
#include "DARWIN_OPTIONS.h"
#endif


CBOP
C !ROUTINE: GCHEM_FORCING_SEP
C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_PREFORMED_TRACERS( dic_offset, bi, bj, 
     &          myTime, myIter, myThid )

C !DESCRIPTION:
C     calls subroutine that will update passive preformed tracer values
C     with a separate timestep. 

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#include "PTRACERS_PARAMS.h"
#include "GCHEM.h"
#ifdef ALLOW_DIC
#include "DIC_VARS.h"
#endif /* ALLOW_DIC */
#ifdef ALLOW_DARWIN
#include "DARWIN_FLUX.h"
#include "DARWIN_SIZE.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter, myThid
      INTEGER dic_offset
      INTEGER bi,bj
CEOP

#ifdef ALLOW_GCHEM
#ifdef GCHEM_SEPARATE_FORCING
#if ((defined GCHEM_PREFORMED) || (defined GCHEM_PREFORMED_MLD))
C!LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER iMin,iMax,jMin,jMax,i,j,k
      INTEGER nTrac
      INTEGER iDIC,iALK,iPO4,iDOP,iDO2,iDFE
      INTEGER iAPRE,iPPRE,iPCO2
      INTEGER kmld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      CHARACTER*(MAX_LEN_FNAM) fldName

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('GCHEM_PREFORMED_TRACERS',myThid)
#endif

      iDIC =UNSET_I
      iALK =UNSET_I
      iPO4 =UNSET_I
      iDOP =UNSET_I
      iDO2 =UNSET_I
      iDFE =UNSET_I
      iAPRE=UNSET_I
      iPPRE=UNSET_I
      iPCO2=UNSET_I
      
      iMin=1-Olx
      iMax=sNx+Olx
      jMin=1-Oly
      jMax=sNy+Oly

C Get arrays of mixed layer depth for preformed tracers
C       DO bj = myByLo(myThid), myByHi(myThid)
C        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=jMin,jMax
          DO i=iMin,iMax
           IF ( usingPCoords ) THEN
            kmld(i,j,bi,bj) = Nr
           ELSE
            kmld(i,j,bi,bj) = 1
           ENDIF
#ifdef GCHEM_PREFORMED_TRACERS_MLD
           DO k=1,Nr
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 
     &           .AND. abs(RF(1+kmld(i,j,bi,bj))) 
     &                .LT.hMixLayer(i,j,bi,bj) 
     &          ) THEN
             kmld(i,j,bi,bj)=kmld(i,j,bi,bj)+1
            ENDIF
           ENDDO
#endif
          ENDDO
         ENDDO
C        ENDDO
C       ENDDO

#ifdef ALLOW_DIC
C get tracer names that are within the dic package
       DO nTrac=1,dic_offset
        fldName=PTRACERS_names(nTrac)
                       
        IF ( fldName(1:3) .EQ. 'dic' ) THEN
            iDIC=nTrac
        ELSEIF ( fldName(1:3) .EQ. 'alk' ) THEN         
            iALK=nTrac
        ELSEIF ( fldName(1:3) .EQ. 'po4' ) THEN  
            iPO4=nTrac   
        ELSEIF ( fldName(1:3) .EQ. 'dop' ) THEN 
            iDOP=nTrac    
        ELSEIF ( fldName(1:2) .EQ. 'o2' ) THEN
            iDO2=nTrac     
        ELSEIF ( fldName(1:2) .EQ. 'fe' ) THEN
            iDFE=nTrac     
        ENDIF
       ENDDO

C JML Do Preformed Tracers associated with the DIC package 
      IF ( PTRACERS_numInUse .GT. dic_offset ) THEN
       DO nTrac=dic_offset+1,PTRACERS_numInUse
        fldName=PTRACERS_names(nTrac)
                       
        IF ( fldName(1:7) .EQ. 'passive' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &      'EVALUATING PASSIVE TRACER',myThid)
#endif 
C Passive tracer set to a constant concentration at the surface (template)        
C            DO bj = myByLo(myThid), myByHi(myThid)
C             DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=jMin,jMax
               DO i=iMin,iMax
                DO k=1,kmld(i,j,bi,bj)
                 pTracer(i,j,k,bi,bj,nTrac)=1. _d 0
                ENDDO
               ENDDO
              ENDDO
C             ENDDO
C            ENDDO
                       
        ELSEIF ( fldName(1:4) .EQ. 'cpre' ) THEN 
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING CPRE TRACER',myThid)
#endif     
C Preformed Carbon, set to surface DIC concentration      
C            DO bj = myByLo(myThid), myByHi(myThid)
C             DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=jMin,jMax
               DO i=iMin,iMax
                DO k=1,kmld(i,j,bi,bj)
                 pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,iDIC)
                ENDDO
               ENDDO
              ENDDO
C             ENDDO
C            ENDDO
            
        ELSEIF ( fldName(1:4) .EQ. 'apre' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING APRE TRACER',myThid)
#endif 
C Preformed Alkalinity, set to surface ALK concentration
C            DO bj = myByLo(myThid), myByHi(myThid)
C             DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=jMin,jMax
               DO i=iMin,iMax
                DO k=1,kmld(i,j,bi,bj)
                 pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,iALK)
                ENDDO
               ENDDO
              ENDDO
C             ENDDO
C            ENDDO
            iAPRE=nTrac 
            
        ELSEIF ( fldName(1:4) .EQ. 'ppre' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING PPRE TRACER',myThid)
#endif 
C Preformed Phosphate, set to surface PO4 concentration    
C            DO bj = myByLo(myThid), myByHi(myThid)
C             DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=jMin,jMax
               DO i=iMin,iMax
                DO k=1,kmld(i,j,bi,bj)
                 pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,iPO4)
                ENDDO
               ENDDO
              ENDDO
C             ENDDO
C            ENDDO
            iPPRE=nTrac   
                      
        ELSEIF ( fldName(1:4) .EQ. 'opre' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING OPRE TRACER',myThid)
#endif 
C Preformed Oxygen, set to surface O2 concentration            
C            DO bj = myByLo(myThid), myByHi(myThid)
C             DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=jMin,jMax
               DO i=iMin,iMax
                DO k=1,kmld(i,j,bi,bj)
                 pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,iDO2)
                ENDDO
               ENDDO
              ENDDO
C             ENDDO
C            ENDDO      
            
        ELSEIF ( fldName(1:6) .EQ.'wm_age') THEN 
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &          'EVALUATING AGE TRACER',myThid)
#endif   
C Watermass Age Tracer, set to zero seconds at the surface           
C            DO bj = myByLo(myThid), myByHi(myThid)
C             DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=jMin,jMax
               DO i=iMin,iMax
                DO k=1,kmld(i,j,bi,bj)
                 pTracer(i,j,k,bi,bj,nTrac)= 0. _d 0
                ENDDO
               ENDDO
              ENDDO
C             ENDDO
C            ENDDO      

        ELSEIF ( fldName(1:7) .EQ. 'atmpco2' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &      'EVALUATING ATMPCO2 TRACER',myThid)
#endif 
C Atmospheric CO2 tracer, set surface to current atmpco2
C   Records what the atmosphere was last time the watermass was at the 
C   surface - may be useful for calculating anthropogenic CO2
C            DO bj = myByLo(myThid), myByHi(myThid)
C             DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=jMin,jMax
               DO i=iMin,iMax
                DO k=1,kmld(i,j,bi,bj)
                 pTracer(i,j,k,bi,bj,nTrac)=AtmospCO2(i,j,bi,bj)*1. _d 6
                ENDDO
               ENDDO
              ENDDO
C             ENDDO
C            ENDDO 
            iPCO2=nTrac
        ENDIF 
       ENDDO
      ENDIF
      
#endif /* ALLOW_DIC */        
#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('GCHEM_PREFORMED_TRACERS',myThid)
#endif
#endif /* GCHEM_PREFORMED */
#endif /* GCHEM_SEPARATE_FORCING */
#endif /* ALLOW_GCHEM */

      RETURN
      END
