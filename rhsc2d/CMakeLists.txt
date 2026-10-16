# 2-D Cartesian geometry, short characteristics

rh_program(rhsc2d rhsc2d
  anglequad.c bezier_2D.c formal.c geometry.c readatmos_xdr.c rhsc2d.c
  piecewise_2D.c piecestokes_2D.c solvelong.c project.c rii.c
  writegeom_xdr.c writeflux_xdr.c)
rh_program(rhsc2d_solve2d solve2d INSTALL_NAME solve2d
  anglequad.c bezier_2D.c iter_2d.c geometry.c piecewise_2D.c
  piecestokes_2D.c solve2d.c solvelong.c project.c)
rh_program(rhsc2d_formal2d formal2d INSTALL_NAME formal2d
  anglequad.c bezier_2D.c formal2d.c geometry.c piecewise_2D.c
  piecestokes_2D.c solvelong.c project.c)
rh_program(rhsc2d_solveray solveray INSTALL_NAME solveray_2d
  anglequad.c bezier_2D.c formal.c geometry.c rii.c readatmos_xdr.c
  solveray.c piecewise_2D.c piecestokes_2D.c solvelong.c project.c)
rh_program(rhsc2d_chartflow chartflow INSTALL_NAME chartflow
  anglequad.c bezier_2D.c chartflow.c formal.c geometry.c readatmos_xdr.c
  rii.c piecewise_2D.c piecestokes_2D.c solvelong.c project.c)
rh_program(rhsc2d_backgrcontr backgrcontr INSTALL_NAME backgrcontr_2d
  anglequad.c backgrcontr.c readatmos_xdr.c project.c)
