!     
!     LOREC  by Xi (Rossi) Luo <xi.rossi.luo@gmail.com>
!     first created on June 13, 2011 
!     Based on previous Matlab version back in 2010
!     
subroutine lorec(n,asig,al,as,dlambda,delta,eps,maxit,ierr)
  !     Fixed gamma, so no line search is performed
  !    delta is a number,  diagonal not penalized
  implicit double precision(a-h,o-z), integer(i-n)
  double precision asig(n,n), al(n,n), as(n,n)
  double precision, dimension(:,:), allocatable :: alold, asold, adev
  parameter(s=2.0)
  ! init arrays
  ierr = 0
  allocate(alold(1:n,1:n), stat=jerr)
  ierr = ierr + jerr
  allocate(asold(1:n,1:n), stat=jerr)
  ierr = ierr + jerr
  allocate(adev(1:n,1:n), stat=jerr)
  ierr = ierr + jerr
  change = 100.0*eps
  v=1.0
  f=0.0
  fold=0.0
  alold = al
  asold = as
  alphaold = 1.0
  alpha = 1.0
  do it=1, maxit
     if (change.lt.eps) exit
     call threshstep(n, asig, alold, asold, adev,al, as, s, dlambda, delta, f)
     if (it.gt.1)  then
        adev = al - alold
        change = frob(n, adev)/(1.0+frob(n, al))
        adev = as - asold
        change = change + frob(n, adev)/(1.0+frob(n, as))
        ! call dblepr("change=",6,change, 1) 
     end if
     fold = f
     alphaold = alpha
     alpha = (1+sqrt(1+4*alphaold*alphaold))*0.5
     v = (alphaold -1)/alpha
     adev = al+v*(al-alold)
     alold = al
     al = adev
     adev = as+v*(as-asold)
     asold = as
     as = adev
  end do
  maxit = it
  if (allocated(adev)) deallocate(adev)
  if (allocated(asold)) deallocate(asold)
  if (allocated(alold)) deallocate(alold)
end subroutine lorec

subroutine threshstep(n, asig, al, as, adev, aln, asn, s, dlambda, delta, f)
  implicit double precision(a-h,o-z), integer(i-n)
  double precision, intent(in) ::  asig(n,n), al(n,n), as(n,n)
  double precision, intent(out) :: adev(n,n), aln(n,n), asn(n,n)
  double precision, dimension(:), allocatable :: ad, awork, aw
  double precision, dimension(:,:), allocatable :: au, atmp
  
  lmax=5*n
  ! init
  allocate(ad(1:n), stat=ierr)
  allocate(aw(1:n), stat=ierr)
  allocate(awork(1:(lmax)), stat=ierr)
  allocate(au(1:n,1:n), stat=ierr)
  allocate(atmp(1:n,1:n), stat=ierr)

  call gradstep(n, asig, al, as, adev, aln, asn, s)
  ! L update
  ! eigen V from lapack  au*aw*au'
  au = aln
  lwork=-1
  call dsyev('V', 'U', n, au, n, aw, awork, lwork, info)
  lwork=min(lmax, int(awork(1)))
  call dsyev('V', 'U', n, au, n, aw, awork, lwork, info)
  idx=0
  sval=0.0
  ! call dblepr("eig=", 4, aw, n)
  do i=1,n
     v=aw(i) - dlambda/s
     if (v.gt.0) then
        atmp(i,:) = v*au(:,i)
        sval = sval + v
        if (idx.eq.0) idx = i
     end if
  end do
  ! call intpr('idx=', 4, idx, 1)
  if (idx.eq.0) then
     aln=0.0
  else
     aln=matmul(au(:,idx:n),atmp(idx:n,:))
  end if
  ! call dblepr("aln12=", 6,  aln(1,2), 1)
  ! call dblepr("aln21=", 6,  aln(2,1), 1)
  
  ! S update
  tmp = delta/s
  dl1 = 0.0
  call softthresh(n, asn, tmp, dl1)

  f = fobj(n, asig, aln, asn) + dlambda*sval + delta*dl1
  
  ! clean up
  if (allocated(atmp)) deallocate(atmp)
  if (allocated(au)) deallocate(au)
  if (allocated(awork)) deallocate(awork)
  if (allocated(aw)) deallocate(aw)
  if (allocated(ad)) deallocate(ad)
end subroutine threshstep



subroutine gradstep(n, asig, al, as, adev, aln, asn, s)
  implicit double precision(a-h,o-z), integer(i-n)
  double precision, intent(in) ::  asig(n,n), al(n,n), as(n,n)
  double precision, intent(out) :: adev(n,n), aln(n,n), asn(n,n)
  call fdev(n, asig, al, as, adev)
  aln=al-(1.0/s)*adev
  asn=as-(1.0/s)*adev
end subroutine gradstep

subroutine fdev(n, asig, al, as, adev)
  implicit double precision(a-h,o-z), integer(i-n)
  double precision, intent(in) ::  asig(n,n), al(n,n), as(n,n)
  double precision, intent(out) :: adev(n,n)
  adev=as+al-asig
end subroutine fdev



subroutine  softthresh(n, x, rho, dl1)
  ! soft threhold matrix x, except diagonal, sum up l1 norm in dl1
  ! assume symmetric x
  implicit double precision(a-h,o-z), integer(i-n)
  double precision, intent(inout) :: x(n, n)
  double precision, intent(in) :: rho 
  dl1 = 0.0
  do i=1,n-1
     do j=(i+1),n
        v=abs(x(i,j))-rho
        if (v.gt.0.0) then 
           x(i,j) = sign(v, x(i,j))
           x(j,i) = x(i,j)
           dl1 = dl1+v*2
        else
           x(i,j) = 0.0
           x(j,i) = 0.0
        end if
     end do
  end do
end subroutine softthresh

double precision  function fobj(n, asig, al, as)
  implicit double precision(a-h,o-z), integer(i-n)
  double precision, intent(in) ::  asig(n,n), al(n,n), as(n,n)
  fobj=0.0
  do i=1,n
     do j=1,n
        fobj = fobj + (asig(i,j)-al(i,j)-as(i,j))**2
     end do
  end do
  fobj=fobj*0.5;
end function fobj

double precision function frob(n, a)
  implicit double precision(a-h,o-z), integer(i-n)
  double precision, intent(in) ::  a(n,n)
  frob =0.0
  do i=1,n
     do j=1,n
        frob = frob + a(i,j)**2
     end do
  end do
  frob = sqrt(frob)
end function frob
