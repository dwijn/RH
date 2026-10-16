# Tools for preparing atomic data. constructatom.c, adjustbf.c and lande.c
# are not built (see issue Q-10).

rh_program(addgrad addgrad addgrad.c modelatom.c rawatom.c)
rh_program(convertatom convertatom convertatom.c modelatom.c)
rh_program(topbase topbase topbase.c modelatom.c)
rh_program(listincl listincl listincl.c)
rh_program(listincl_f listincl_f listincl_f.c)
rh_program(printneff printneff printneff.c rawatom.c)
rh_program(impact impact impact.c rawatom.c laguerre.c bessel.c)
rh_program(hcollisions hcollisions hcollisions.c laguerre.c)
rh_program(make_h make_h make_h.c modelatom.c)
rh_program(avgmolweight avgmolweight avgmolweight.c)
rh_program(readnist readnist readnist.c)
rh_program(nist_pf nist_pf nist_pf.c)
rh_program(waveinfo waveinfo waveinfo.c)
rh_program(readkurucz readkurucz readkurucz.c writemulti.c)
rh_program(readtrujillo readtrujillo readtrujillo.c writemulti.c)
