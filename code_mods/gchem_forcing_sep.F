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
      SUBROUTINE GCHEM_FORCING_SEP( myTime, myIter, myThid )

C !DESCRIPTION:
C     calls subroutine that will update passive tracers values
C     with a separate timestep. Since GCHEM_FORCING_SEP is now
C     called before DO_FIELDS_BLOCKING_EXCHANGES, the passive
C     tracer values in the halo regions are not up to date and
C     must not be used.

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
CEOP

#ifdef ALLOW_GCHEM
#ifdef GCHEM_SEPARATE_FORCING

C!LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER bi,bj,iMin,iMax,jMin,jMax
c     INTEGER i,j,
      INTEGER dic_offset

      dic_offset=UNSET_I
      
ccccccccccccccccccccccccc
c global calculations   c
ccccccccccccccccccccccccc
#ifdef ALLOW_OLD_VIRTUALFLUX
#ifdef ALLOW_DIC
c find global surface averages
       gsm_s = 0. _d 0
       gsm_dic = 0. _d 0
       gsm_alk = 0. _d 0
       CALL TRACER_MEANAREA(salt, 1,gsm_s,myThid)
       CALL TRACER_MEANAREA(
     &             pTracer(1-OLx,1-OLy,1,1,1,1), 1, gsm_dic, myThid )
       print*,'mean surface dic', gsm_dic,gsm_s
       CALL TRACER_MEANAREA(
     &             pTracer(1-OLx,1-OLy,1,1,1,2), 1, gsm_alk, myThid )
#endif
#ifdef ALLOW_DARWIN
c find global surface averages
       gsm_s = 0. _d 0
       gsm_dic = 0. _d 0
       gsm_alk = 0. _d 0
       CALL TRACER_MEANAREA(salt, 1,gsm_s,myThid)
       CALL TRACER_MEANAREA(
     &             pTracer(1-OLx,1-OLy,1,1,1,iDIC), 1, gsm_dic, myThid )
       print*,'mean surface dic', gsm_dic,gsm_s
       CALL TRACER_MEANAREA(
     &             pTracer(1-OLx,1-OLy,1,1,1,iALK), 1, gsm_alk, myThid )
#endif
#endif /* ALLOW_OLD_VIRTUALFLUX */
ccccccccccccccccccccccccccccccccccccccccccc

#ifdef ALLOW_DARWIN
       IF ( useDARWIN ) THEN
         CALL DARWIN_CONS( myIter, myTime, myThid )
       ENDIF
#endif

ccccccccccccccccccccccccc
c chemical forcing      c
ccccccccccccccccccccccccc
C$taf loop = parallel
       DO bj=myByLo(myThid),myByHi(myThid)
C$taf loop = parallel
        DO bi=myBxLo(myThid),myBxHi(myThid)

        jMin=1
        jMax=sNy
        iMin=1
        iMax=sNx

ccccccccccccccccccccccccccc DIC cccccccccccccccccccccccccccccccc

#ifdef ALLOW_DIC
#ifdef ALLOW_FE
#ifdef ALLOW_VARIABLE_LIGANDS
          dic_offset=7
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,6),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,7),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
#else
          dic_offset=6
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,6),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
#endif /* ALLOW_VARIABLE_LIGANDS */
#else
#ifdef ALLOW_O2
          dic_offset=5
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
#else
          dic_offset=4
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
#endif /* ALLOW_O2 */
#endif /* ALLOW_FE */
#endif /* ALLOW_DIC */

cccccccccccccccccccccccccc END DIC cccccccccccccccccccccccccccccccccc
#ifdef ALLOW_DARWIN
        IF ( useDARWIN ) THEN
#ifdef NUT_SUPPLY
c articficial supply of nutrients
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
#endif
ccccccccccccccc
C darwin_forcing operates on bi,bj part only, but needs to get full
C array because of last (iPtr) index
          CALL DARWIN_FORCING(  pTracer(1-OLx,1-OLy,1,1,1,1),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
C JML dont know if this is even right?     
          dic_offset=ngud
        ENDIF
#endif /* ALLOW_DARWIN */

        ENDDO
       ENDDO

#ifdef ALLOW_DARWIN
       IF ( useDARWIN ) THEN
         CALL DARWIN_CONS( myIter, myTime, myThid )
#ifdef ALLOW_CARBON
         CALL DIC_ATMOS( 1, myTime, myIter, myThid )
#endif
       ENDIF
#endif /* ALLOW_DARWIN */

#ifdef ALLOW_DIC
       CALL DIC_ATMOS( 1, myTime, myIter, myThid )
       CALL DIC_STORE_FLUXCO2( myTime, myIter, myThid )
#endif

#if ((defined GCHEM_PREFORMED) || (defined GCHEM_PREFORMED_MLD))
C JML If the number of ptracers is greater than that required by 
C  pkg DIC (or DARWIN) then run this routine for the preformed tracers.
       IF (PTRACERS_numInUse.GT.dic_offset) THEN
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          CALL GCHEM_PREFORMED_TRACERS( dic_offset, bi, bj,  
     &                    myTime, myIter, myThid )
         ENDDO
        ENDDO
       ENDIF
#endif

#ifdef ALLOW_COST
       CALL DIC_COST( myTime, myIter, myThid )
#endif

#endif /* GCHEM_SEPARATE_FORCING */
#endif /* ALLOW_GCHEM */

      RETURN
      END
