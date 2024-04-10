#define ESMF_ERR_ABORT(rc) if (ESMF_LogFoundError(rc, msg="Aborting UFS", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
program convertics

  ! TODO: r8 should be specified

  use ESMF

  use utils_mod, only : getfield, getvecpair, dumpnc, remap, nf90_err, chkerr
  use init_mod
  use netcdf

  implicit none

  integer, parameter :: srcnx = 1440, srcny = 1080, ncat = 5, nilyr = 7
  integer, parameter :: dstnx = 360 , dstny =  320

  character(len=140) :: dirgrd = '/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data/'
  character(len=140) :: dirsrc = '/scratch1/NCEPDEV/stmp2/Denise.Worthen/tripole/fix.mx025/RESTART/'
  character(len=140) :: fwgt, input_file

  real(kind=8), dimension(srcnx*srcny) :: angsrc, cosrotsrc, sinrotsrc
  real(kind=8), dimension(dstnx*dstny) :: angdst, cosrotdst, sinrotdst

  real(kind=8), dimension(srcnx*srcny,2) :: vecpairsrc, tmpsrc
  real(kind=8), dimension(dstnx*dstny,2) :: vecpairdst, tmpdst

  integer :: ii, ncid, rc
  real(kind=8) :: urot, vrot

  type(ESMF_VM) :: vm
  type(ESMF_RouteHandle) :: rh
  type(ESMF_Mesh)  :: meshsrc, meshdst
  type(ESMF_Field) :: fldsrc, flddst

  real(kind=8) , pointer  :: srcptr(:), dstptr(:)

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

  call readnml

  ! get the rotation angles on Ct
  input_file = trim(dirgrd)//trim(fsrc)//'/tripole.mx'//trim(fsrc)//'.nc'
  print *,trim(input_file)
  call nf90_err(nf90_open(trim(input_file), nf90_nowrite, ncid), 'open: '//trim(input_file))
  call getfield(trim(input_file), 'anglet', dims=(/srcnx,srcny/), field=angsrc)
  call nf90_err(nf90_close(ncid), 'close: '//trim(input_file))
  ! reverse direction for CICE
  cosrotsrc =  cos(-angsrc)
  sinrotsrc =  sin(-angsrc)

  input_file = trim(dirgrd)//trim(fdst)//'/tripole.mx'//trim(fdst)//'.nc'
  print *,trim(input_file)
  call nf90_err(nf90_open(trim(input_file), nf90_nowrite, ncid), 'open: '//trim(input_file))
  call getfield(trim(input_file), 'anglet', dims=(/dstnx,dstny/), field=angdst)
  call nf90_err(nf90_close(ncid), 'close: '//trim(input_file))
  ! reverse direction for CICE
  cosrotdst =  cos(-angdst)
  sinrotdst =  sin(-angdst)

  ! get the src velocities, on Ct and rotate EW
  input_file = trim(dirsrc)//'iced.2011-10-03-00000.nc'
  print *,trim(input_file)
  call getvecpair(trim(input_file), trim(dirgrd)//trim(fsrc)//'/', cosrotsrc, sinrotsrc, 'uvel', 'Bu', &
       'vvel', 'Bu', dims=(/srcnx,srcny/), vecpair=vecpairsrc)
  print *,' got vecpair'

  !print *,vecpairsrc(:,1),vecpairsrc(:,2)
  !map to Ct
  fwgt = trim(dirgrd)//trim(fsrc)//'/tripole.mx'//trim(fsrc)//'.Bu.to.Ct.bilinear.nc'
  !tmpsrc = 0.0
  !call remap(trim(fwgt), dim2=2, src_field=vecpairdst, dst_field=tmpsrc)
  !vecpairsrc = tmpsrc
  call dumpnc('testsrc.Ct_ew.nc','vecsrc', dims=(/srcnx,srcny/), nflds=2, field=vecpairsrc)
  !------
  !------

  call ESMF_Initialize(rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_VMGetGlobal(vm, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  input_file = trim(dirgrd)//trim(fsrc)//'/mesh.mx'//trim(fsrc)//'.nc'
  print *,trim(input_file)
  meshsrc = ESMF_MeshCreate(filename=trim(input_file), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  fldsrc = ESMF_FieldCreate(meshsrc, ESMF_TYPEKIND_R8, name='fsrc', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  input_file = trim(dirgrd)//trim(fdst)//'/mesh.mx'//trim(fdst)//'.nc'
  print *,trim(input_file)
  meshdst = ESMF_MeshCreate(filename=trim(input_file), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  flddst = ESMF_FieldCreate(meshdst, ESMF_TYPEKIND_R8, name='fdst', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  print *,'call regrid'

  call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=rh, &
       srcMaskValues=(/0/),                                  &
       dstMaskValues=(/0/),                                  &
       regridmethod=ESMF_REGRIDMETHOD_BILINEAR,              &
       extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD,          &
       polemethod=ESMF_POLEMETHOD_ALLAVG,                    &
       !srcTermProcessing=srcTermProcessing_Value,           &
       ignoreDegenerate=.true.,                              &
       !dstStatusField=dststatusfield,                       &
       unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! get pointers to src and dst fields and fill src side
  call ESMF_FieldGet(fldsrc, farrayPtr=srcptr, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_FieldGet(flddst, farrayPtr=dstptr, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! remap uvel
  dstptr(:) = 0.0
  srcptr(:) = vecpairsrc(:,1)
  call ESMF_FieldRegrid(fldsrc, flddst, routehandle=rh, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  vecpairdst(:,1) = dstptr(:)

  ! remap vvel
  dstptr(:) = 0.0
  srcptr(:) = vecpairsrc(:,2)
  call ESMF_FieldRegrid(fldsrc, flddst, routehandle=rh, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  vecpairdst(:,2) = dstptr(:)
  !call dumpnc('testdstCt_ew.nc','vecdst', dims=(/dstnx,dstny/), nflds=2, field=vecpairdst)

  ! rotate back to ij on Ct
  urot = 0.0; vrot = 0.0
  do ii = 1,dstnx*dstny
     urot = vecpairdst(ii,1)*cosrotdst(ii) + vecpairdst(ii,2)*sinrotdst(ii)
     vrot = vecpairdst(ii,2)*cosrotdst(ii) - vecpairdst(ii,1)*sinrotdst(ii)
     vecpairdst(ii,1) = urot
     vecpairdst(ii,2) = vrot
  end do
  call dumpnc('testdstCt_ij.nc','vecdst', dims=(/dstnx,dstny/), nflds=2, field=vecpairdst)

  ! stagger back to Bu
  fwgt = trim(dirgrd)//trim(fdst)//'/tripole.mx'//trim(fdst)//'.Ct.to.Bu.bilinear.nc'
  print *,trim(fwgt)
  tmpdst = 0.0
  call remap(trim(fwgt), dim2=2, src_field=vecpairdst, dst_field=tmpdst)
  print *,'XXX ',minval(tmpdst),maxval(tmpdst)
  call dumpnc('testdstBu_ij.nc','vecdst', dims=(/dstnx,dstny/), nflds=2, field=tmpdst)

end program convertics
