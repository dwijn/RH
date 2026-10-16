# 1-D plane-parallel geometry

rh_program(rhf1d rhf1d
  anglequad.c bezier_1D.c feautrier.c formal.c hydrostat.c multiatmos.c
  rhf1d.c piecewise_1D.c project.c riiplane.c writeflux_xdr.c
  writegeom_xdr.c)
rh_program(rhf1d_solve1d solve1d INSTALL_NAME solve1d
  bezier_1D.c solve1d.c iter_1d.c feautrier.c piecewise_1D.c project.c)
rh_program(rhf1d_solveray solveray INSTALL_NAME solveray_1d
  anglequad.c bezier_1D.c feautrier.c formal.c hydrostat.c multiatmos.c
  project.c piecewise_1D.c riiplane.c solveray.c writegeom_xdr.c)
rh_program(rhf1d_backgrcontr backgrcontr INSTALL_NAME backgrcontr_1d
  anglequad.c backgrcontr.c multiatmos.c project.c)
