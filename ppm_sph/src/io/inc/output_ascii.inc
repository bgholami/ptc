  do i = prtl_1, prtl_N

    !< Plot only desired phases
    if (phase(ap(phaseidx,i))%output .eq. 0) cycle

    rfac = 1._mk
#ifdef __NONDIMOUTPUT
    rfac = 1._mk/length_ref
#endif
    do ii = 1, ndim
      outvec(ii) = xp(ii,i)*rfac
    enddo
    ldaoutvec_ascii = ndim

    !< Velocity
    if (output_velocity) then
      rfac = 1._mk
#ifdef __NONDIMOUTPUT
      rfac = 1._mk/vel_ref
#endif
      do ii = 1, ndim
        outvec(ldaoutvec_ascii+ii) = dxpdt(ii,i)*rfac
      enddo
      ldaoutvec_ascii = ldaoutvec_ascii + ndim
    endif

    !< Density
    if (output_density) then
      rfac = 1._mk
#ifdef __NONDIMOUTPUT
      rfac = 1._mk/density_ref
#endif
      outvec(ldaoutvec_ascii+1) = pdata(rhoidx,i)*rfac
      ldaoutvec_ascii = ldaoutvec_ascii + 1
    endif

    !< Pressure
    if (output_pressure) then
      rfac = 1._mk
#ifdef __NONDIMOUTPUT
      rfac = 1._mk/(density_ref*vel_ref**2)
#endif
      outvec(ldaoutvec_ascii+1) = pdata(pressidx,i)*rfac
      ldaoutvec_ascii = ldaoutvec_ascii + 1
    endif

    !< Volume
    if (output_volume) then
      rfac = 1._mk
#ifdef __NONDIMOUTPUT
      rfac = 1._mk/(length_ref**__DIM)
#endif
      outvec(ldaoutvec_ascii+1) = pdata(volidx,i)*rfac
      ldaoutvec_ascii = ldaoutvec_ascii + 1
    endif

    !< Volume
    if (output_mass) then
      rfac = 1._mk
#ifdef __NONDIMOUTPUT
      rfac = 1._mk/(density_ref*length_ref**__DIM)
#endif
      outvec(ldaoutvec_ascii+1) = pdata(massidx,i)*rfac
      ldaoutvec_ascii = ldaoutvec_ascii + 1
    endif

    if (output_h) then
      rfac = 1._mk
#ifdef __NONDIMOUTPUT
      rfac = 1._mk/length_ref
#endif
      outvec(ldaoutvec_ascii+1) = pdata(hidx,i)*rfac
      ldaoutvec_ascii = ldaoutvec_ascii + 1
    endif

    if (output_phaseidx) then
      outvec(ldaoutvec_ascii+1) = real(ap(phaseidx,i),mk)
      ldaoutvec_ascii = ldaoutvec_ascii + 1
    endif
    if (output_id) then
      outvec(ldaoutvec_ascii+1) = real(ap(partIDidx,i),mk)
      ldaoutvec_ascii = ldaoutvec_ascii + 1
    endif
    if (output_rank) then
      outvec(ldaoutvec_ascii+1) = real(rank,mk)
      ldaoutvec_ascii = ldaoutvec_ascii + 1
    endif

    call ppm_io(IUNIT=unit_writeascii,ADATA=outvec,ACTN=ppm_param_io_write,&
    DIST=ppm_param_io_concat,IOFMT=trim(iofmt),STAT=info)
    if (info .ne. 0) call abortmessage(myname,'Writing '//trim(filename)//' failed..')
  enddo
