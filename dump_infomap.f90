! Subroutine of PATHSAMPLE to dump rate constant and equilibrium probability information for use with
! e.g. InfoMap or MLR-MCL algorithms
! djs244

SUBROUTINE DUMP_INFOMAP

    USE COMMONS, ONLY: PFMIN, EMIN, TEMPERATURE, FVIBMIN, HORDERMIN, &
                       KPLUS, KMINUS, NMIN, NTS

    IMPLICIT NONE

    DOUBLE PRECISION :: PEQ(NMIN), PEQ_TOT, PFNORM, PFMEAN
    INTEGER :: J1

    PRINT *, 'dump_infomap> dumping rate constant information'

    OPEN(6,FILE="ts_weights.dat",STATUS="new")
    ! logarithms of rate constants: 2 lines (x->y and y->x) for each line in ts.data
    DO J1=1,NTS
        WRITE(6,FMT="(F64.58)") KPLUS(J1)
        WRITE(6,FMT="(F64.58)") KMINUS(J1)
    ENDDO
    CLOSE(6)

    PRINT *, 'dump_infomap> dumping equilibrium probabilities for minima'

    PFNORM=0.0D0
    PFMEAN=0.0D0
    DO J1=1,NMIN
        PFMIN(J1) = -EMIN(J1)/TEMPERATURE - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
!        PFNORM=PFNORM+EXP(PFMIN(J1))
        PFMEAN = PFMEAN+PFMIN(J1) 
    ENDDO
!    PFNORM=LOG(PFNORM)
    PFMEAN = PFMEAN/FLOAT(NMIN)

    DO J1=1,NMIN
        PFMIN(J1) = -EMIN(J1)/TEMPERATURE - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
        PFNORM=PFNORM+EXP(PFMIN(J1)-PFMEAN)
    ENDDO
    PFNORM = LOG(PFNORM)+PFMEAN

    PEQ_TOT = 0.D0
    DO J1=1,NMIN
        PEQ(J1)=EXP(PFMIN(J1)-PFNORM)
        PEQ_TOT = PEQ_TOT + PEQ(J1)
    ENDDO
    PRINT '(A,G20.10)','sum of equilibrium occupation probabilities=',PEQ_TOT

    ! print logarithms of equilibrium occupation probabilities of minima
    OPEN(7,FILE="stat_prob.dat",STATUS="new")
    DO J1=1,NMIN
        WRITE(7,FMT="(F64.58)") LOG(PEQ(J1))
    ENDDO
    CLOSE(7)

END SUBROUTINE DUMP_INFOMAP
