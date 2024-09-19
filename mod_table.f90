Module mod_table

implicit none

! ========================================================================

  character(2)                    :: Table(103)
  character(5)                    :: calc_flags(10)
  character(1)                    :: Orbital(0:6)
  integer                         :: OrbSize(0:6)
  character(10)                   :: labelSph(49)
  integer                         :: labSphOff(-1:6)

  data Table /' H',                                                                                'He',  &
              'Li','Be',                                                  ' B',' C',' N',' O',' F','Ne',  &
              'Na','Mg',                                                  'Al','Si',' P',' S','Cl','Ar',  &
              ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',  &
              'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe',  &
              'Cs','Ba',                                                                                  &
                        'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',       &
                             'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',  &
              'Fr','Ra',                                                                                  &
                        'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr' /

  data Orbital /'S','P','D','F','G','H','I'/

  data OrbSize /1, 3, 5, 7, 9, 11, 13/

  data labSphOff /0, 1,  4,  9, 16, 25, 36, 49/

  data labelSph /                                                             '         s',                                                                              &
                                                                 '        px','        py','        pz',                                                                 &
                                                    '       dxy','       dyz','       dz2','       dxz','    dx2-y2',                                                    &
                                       'fy(3x2-y2)','      fxyz','      fyz2','       fz3','      fxz2','      fzx2','       fx3',                                       &
                          '      gyx3','     gx2yz','     gxyz2','      gyz3','       gz4','      gxz3','     gx2z2','      gzx3','       gy4',                          &
             '      hyx4','     hx3yz','    hyx2z2','     hxyz3','      hyz4','    hx2y2z','      hxz4','     hx2z3','     hx3z2','      hzy4','      hxy4',             &
'       i-6','       i-5','       i-4','       i-3','       i-2','       i-1','        i0','        i1','        i2','        i3','        i4','        i5','        i6' /

  data calc_flags /'dimer','trime','sp   ','opt  ','scan ','genet','full ','simpl','mopri','nomol'/ 
  !Write the flags in LOWERCASE
  !DIMER TRIMEr Single-Point OPT SCAN FULL SIMPLe MOPRInt NOMOLden
  !May shorten them after the program is done

End Module mod_table
