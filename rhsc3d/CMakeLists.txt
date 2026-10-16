# 3-D Cartesian geometry, short characteristics

rh_program(rhsc3d rhsc3d
  anglequad.c bezier_3D.c formal.c geometry.c readatmos_xdr.c rhsc3d.c
  piecewise_3D.c piecestokes_3D.c solvelong.c project.c rii.c
  writeflux_xdr.c writegeom_xdr.c)
rh_program(rhsc3d_solveray solveray INSTALL_NAME solveray_3d
  anglequad.c bezier_3D.c formal.c geometry.c readatmos_xdr.c rii.c
  piecewise_3D.c piecestokes_3D.c solvelong.c solveray.c project.c)
rh_program(rhsc3d_formal3d formal3d INSTALL_NAME formal3d
  anglequad.c bezier_3D.c formal3d.c geometry.c piecewise_3D.c
  piecestokes_3D.c solvelong.c project.c)
