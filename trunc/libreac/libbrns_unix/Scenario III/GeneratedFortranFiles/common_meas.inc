        parameter (maxxmeas=1,maxspmeas=1)
        common/meas1/idspmeas,nrxmeas,nrspmeas
        common/meas2/spmeas,xmeas,sigmeas
        dimension idspmeas(ntopt,maxspmeas),nrxmeas(ntopt,maxxmeas),nrsp
     +meas(ntopt),spmeas(ntopt,maxspmeas,maxxmeas),xmeas(ntopt,maxspmeas
     +,maxxmeas),sigmeas(ntopt,maxspmeas,maxxmeas)
        integer idspmeas,nrxmeas,nrspmeas
        real*8 spmeas,xmeas,sigmeas
