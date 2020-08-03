#include "DIC_OPTIONS.h"
CBOP
C     !ROUTINE: EXPINT
C     !INTERFACE:
      _RL FUNCTION EXPINT(n,x,myThid)
C     !DESCRIPTION: \bv
C     *=============================================================*
C     | S/R  EXPINT(n,x)
C     | o compute exponential integral for "gamma" function remin profile
C     | o evaluates the exponential integral En(x).
C     *=============================================================*
C
C     REFERENCES:
C     NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING 
C             (ISBN 0-521-43064-X) 6.3 Exponential Integrals p 217 
C     CAEL AND BISSON (2018), FMARS, 5, 395, doi:10.3389/fmars.2018.00395
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DIC_VARS.h"
C     !INPUT/OUTPUT PARAMETERS:
C     N is the 
C     MAXIT is the maximum allowed number of iterations; 
C     EPS is the desired relative error, not smaller than the machine precision; 
C     FPMIN is a number near the smallest representable floating-point number;
C     EULER is Euler’s constant.
      INTEGER n, myThid
      _RL x
#if (defined ALLOW_PTRACERS && defined DIC_BIOTIC)      
      INTEGER MAXIT
      _RL EPS,FPMIN,EULER
      PARAMETER (MAXIT=100)
      PARAMETER (EPS=1. _d -7)
      PARAMETER (FPMIN=1. _d -30)
      PARAMETER (EULER=.5772156649 _d 0)

C Parameters: 
       INTEGER i,ii,nm1
       _RL a,b,c,d,del,fact,h,psi
       CHARACTER*(MAX_LEN_MBUF) msgBuf

C Initialize result
       expint=0. _d 0

       nm1=n-1
       
       IF(n.lt.0.or.x.lt.0. _d 0.or.
     &     (x.eq.0. _d 0.and.(n.eq.0.or.n.eq.1))) THEN
C            pause ’bad arguments in expint’
             WRITE(msgBuf,'(A)')
     &          'PHOS_FLUX: bad arguments in EXPINT'
             CALL PRINT_ERROR( msgBuf, myThid )
C            STOP 'ABNORMAL END: PHOS_FLUX error' 
       ELSEIF(n.eq.0) THEN 
C      Special case.
           expint=exp(-x)/x
       ELSEIF(x.eq.0.) THEN
C      Another special case.
           expint=1./nm1
       ELSEIF(x.gt.1.) THEN
C      Lentz’s algorithm (§5.2).
           b=x+n
           c=1./FPMIN
           d=1./b
           h=d
           DO i=1,MAXIT
               a=-i*(nm1+i)
               b=b+2.
               d=1./(a*d+b) 
C      Denominators cannot be zero.
               c=b+a/c
               del=c*d
               h=h*del
       
               IF (abs(del-1.).lt.EPS) THEN
                   expint=h*exp(-x)
                   RETURN
               ENDIF
           ENDDO
C           pause ’continued fraction failed in expint’
           WRITE(msgBuf,'(A)')
     &          'PHOS_FLUX: continued fraction failed in EXPINT'
           CALL PRINT_ERROR( msgBuf, myThid )
C           STOP 'ABNORMAL END: PHOS_FLUX error' 
       ELSE
C      Evaluate series.
           IF (nm1.ne.0) THEN
C      Set first term.
               expint=1./nm1
           ELSE
               expint=-log(x)-EULER
           ENDIF
           fact=1.
           DO i=1,MAXIT
               fact=-fact*x/i
               IF (i.ne.nm1) THEN
                   del=-fact/(i-nm1)
               ELSE
                   psi=-EULER 
C                   Compute ψ(n).
                   DO ii=1,nm1
                       psi=psi+1./ii
                   ENDDO
                   del=fact*(-log(x)+psi)
               ENDIF
               expint=expint+del
               IF (abs(del).lt.abs(expint)*EPS) return
           ENDDO
C       pause ’series failed in expint’
           WRITE(msgBuf,'(A)')
     &          'PHOS_FLUX: series failed in EXPINT'
           CALL PRINT_ERROR( msgBuf, myThid )
C           STOP 'ABNORMAL END: PHOS_FLUX error' 
       ENDIF
#endif
       RETURN
       END
