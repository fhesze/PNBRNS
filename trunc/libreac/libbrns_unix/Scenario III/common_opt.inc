        parameter (nopt=1,ntopt=1,ntotparam=1)
        parameter (idoopt=0)
        common/copt1/par,timemeas
        common/copt2/idpar
        dimension par(ntotparam),idpar(nopt),timemeas(ntopt)
        real*8 par,timemeas
        integer idpar
