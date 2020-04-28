

! ====================================================== !
! === linear Interpolation for 1D profile            === !
! ====================================================== !

subroutine LinearInterp1D( xa, fa, xp, fp, nData, nIntp, force_to_ )
  implicit none
  integer         , intent(in)  :: nData, nIntp
  double precision, intent(in)  :: xa(nData), fa(nData)
  double precision, intent(in)  :: xp(nData)
  double precision, intent(out) :: fp(nData)
  integer         , intent(in)  :: force_to_
  integer                       :: ik, is, iIni, iFin, iMid
  double precision              :: dxInv, p1, p2, hxp

  iMid = int( nData / 2 )
  do ik=1, nIntp

     hxp = xp(ik)

     ! ------------------------------------------------------ !
     ! --- [1] boundary value exceptions                  --- !
     ! ------------------------------------------------------ !
     !  -- out-of-range exception :: xp < min(xa)         --  !
     if ( hxp.lt.xa(    1) ) then
        if ( force_to_.eq.1 ) then
           hxp = xa(1)
        else
           write(6,*) "xp is out of range :: ik = ", ik, " xp(ik) = ", hxp
           write(6,*) "             range :: xa(1) = ", xa(1), ", xa(nData) = ", xa(nData)
           stop
        endif
     endif
     !  -- out-of-range exception :: xp > max(xa)         --  !
     if ( hxp.gt.xa(nData) ) then
        if ( force_to_.eq.1 ) then
           hxp = xa(nData)
        else
           write(6,*) "xp is out of range :: ik = ", ik, " xp(ik) = ", hxp
           write(6,*) "             range :: xa(1) = ", xa(1), ", xa(nData) = ", xa(nData)
           stop
        endif
     endif
     !  -- exact boundary value case ::                   --  !
     if ( hxp.eq.xa(nData) ) then
        fp(ik) = fa(nData)
        exit
     endif
     ! ------------------------------------------------------ !
     ! --- [2] search section to be interpolated          --- !
     ! ------------------------------------------------------ !
     !  --  initiarl search position                      --  !
     if ( hxp.ge.xa(iMid) ) then
        iIni = iMid
        iFin = nData
     else
        iIni = 1
        iFin = iMid-1
     endif
     !  --  search & interpolate                          --  !
     do is=iIni, iFin
        if ( ( hxp.ge.xa(is) ).and.( hxp.lt.xa(is+1) ) ) then
           dxInv  = 1.d0 / ( xa(is+1) - xa(is) )
           p1     = ( xa(is+1) - hxp    ) * dxInv
           p2     = ( hxp      - xa(is) ) * dxInv
           fp(ik) = p1*fa(is) + p2*fa(is+1)
           exit
        endif
     enddo

  enddo

  return
end subroutine LinearInterp1D


! ====================================================== !
! === biLinear  Interpolation                        === !
! ====================================================== !
subroutine LinearInterp2D( gData, pData, xMin, yMin, dx, dy, LI, LJ, nData  )
  implicit none
  integer         , intent(in)  :: LI, LJ, nData
  double precision, intent(in)  :: dx, dy, xMin, yMin
  double precision, intent(in)  :: gData(LI,LJ)
  double precision, intent(out) :: pData(3,nData)
  integer                       :: ip, jp, ii, jj, ir, jr, iData
  double precision              :: xp, yp, dxInv, dyInv, sum, xMax, yMax
  double precision              :: sfx(-2:2), sfy(-2:2), sfz(-2:2)
  integer         , parameter   :: x_=1, y_=2, v_=3

  xMax  = dx*dble(LI-1) + xMin
  yMax  = dy*dble(LJ-1) + yMin
  dxInv = 1.d0 / dx
  dyInv = 1.d0 / dy
  do iData=1, nData
     xp   = min( max( pData(x_,iData)-xMin, 0.d0 ), xMax-xMin )
     yp   = min( max( pData(y_,iData)-yMin, 0.d0 ), yMax-yMin )
     ip   = nint( xp*dxInv )
     jp   = nint( yp*dyInv )
     sfx  = shapef( xp, ip, dxInv )
     sfy  = shapef( yp, jp, dyInv )
     sum  = 0.d0
     do jj=-2, 2
        do ii=-2, 2
           ir  = ip+1 + ii
           jr  = jp+1 + jj
           ir  = min( max( 1, ir ), LI )
           jr  = min( max( 1, jr ), LJ )
           sum = sum + sfx(ii)*sfy(jj) * gData(ir,jr)
        enddo
     enddo
     pData(v_,iData) = sum
  enddo

  return
contains
  
  ! ====================================================== !
  ! === shape Function                                 === !
  ! ====================================================== !
  function shapef( rp, lp, drInv )
    implicit none
    integer         , intent(in) :: lp
    double precision, intent(in) :: rp, drInv
    double precision             :: delta
    double precision             :: shapef(-2:+2)

    delta     = rp*drInv - dble(lp)
    shapef(-2) = 0.d0
    shapef(-1) = 0.d0
    shapef( 0) = 1.d0 - abs( delta )
    shapef(+1) = 0.d0
    shapef(+2) = 0.d0
    shapef( int( sign( 1.d0, delta ) ) ) = abs( delta )

    return
  end function shapef

end subroutine LinearInterp2D


