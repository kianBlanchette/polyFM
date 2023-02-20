subroutine modF(mxdeg, nord, FPcoef, Fcoef, Fp, nfreqs)
  implicit none
  real*8 FPcoef(mxdeg,nfreqs), Fcoef(nord,nfreqs) &
    , Fp(mxdeg,nfreqs)
  integer mxdeg, nord, nfreqs, i

  do i=1, nfreqs
    call modFcomp(mxdeg,nord,FPcoef(:,i),Fcoef(:,i),Fp(:,i))
  end do
end subroutine
