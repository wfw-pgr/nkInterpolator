! ====================================================== !
! === linear interpolation form a triangle           === !
! ====================================================== !
subroutine interpolate__fromTriMesh( nodes, elems, points, nNodes, nElems, nPoints )
  implicit none
  integer         , intent(in)    :: nNodes, nElems, nPoints
  integer         , intent(in)    :: elems(3,nElems)
  double precision, intent(in)    :: nodes(3,nNodes)
  double precision, intent(inout) :: points(3,nPoints)
  integer                         :: iN, iP
  double precision                :: vertex(3,3)
  integer         , allocatable   :: enclose(:,:)
  double precision, allocatable   :: node2D(:,:)
  integer         , parameter     :: x_=1, y_=2, z_=3
  integer         , parameter     :: nd1_=1, nd2_=2, nd3_=3
  integer         , parameter     :: elm_=1, edg_=2
  character(20)                   :: interpolation__type = "triangle"

  ! ------------------------------------------------------ !
  ! --- [1] initialization                             --- !
  ! ------------------------------------------------------ !
  allocate( node2D(3,nNodes), enclose(2,nPoints) )
  do iN=1, nNodes
     node2D(x_,iN)  = nodes(x_,iN)
     node2D(y_,iN)  = nodes(y_,iN)
     node2D(z_,iN)  = 0.d0
  enddo
  do iP=1, nPoints
     points(z_,iP) = 0.d0
  enddo
  call investigate__inside_triangle( node2D, elems, points, enclose, nNodes, nElems, nPoints )
  if ( trim(interpolation__type).eq."barycentric" ) then
     do iP=1, nPoints
        if ( enclose(elm_,iP) > 0 ) then
           call barycentric__interpolate( nodes(x_:z_,elems(nd1_,enclose(elm_,iP))), &
                &                         nodes(x_:z_,elems(nd2_,enclose(elm_,iP))), &
                &                         nodes(x_:z_,elems(nd3_,enclose(elm_,iP))), &
                &                         points(:,iP) )
        else
           do iN=nd1_,nd3_
              vertex(x_:z_,iN) = nodes( x_:z_,elems(iN,abs(enclose(elm_,iP)) ) )
           enddo
           call interpolate__trianglePoint( vertex, points(x_:z_,iP) )
        endif
     enddo
  endif

  if ( trim(interpolation__type).eq."triangle" ) then
     do iP=1, nPoints
        do iN=nd1_,nd3_
           vertex(x_:z_,iN) = nodes( x_:z_,elems(iN,abs(enclose(elm_,iP)) ) )
        enddo
        call interpolate__trianglePoint( vertex, points(x_:z_,iP) )
     enddo
  end if

  return
end subroutine interpolate__fromTriMesh


! ====================================================== !
! === investigate__inside_triangle                   === !
! ====================================================== !
subroutine investigate__inside_triangle( nodes, elems, points, enclose, nNodes, nElems, nPoints )
  implicit none
  integer         , intent(in)    :: nNodes, nElems, nPoints
  integer         , intent(in)    :: elems(3,nElems)
  double precision, intent(in)    :: nodes(3,nNodes), points(3,nPoints)
  integer         , intent(inout) :: enclose(2,nPoints)
  integer                         :: iP, iE, iv
  double precision                :: norm_length, inner(3), dist
  double precision                :: normc(3), cog(3)
  double precision                :: edge(3,3), VtoP(3,3), VtoG(3,3), crossP(3,3), crossG(3,3)
  double precision, allocatable   :: distance(:)
  integer         , parameter     :: x_  =1, y_  =2, z_  =3
  integer         , parameter     :: e1_ =1, e2_ =2, e3_ =3
  integer         , parameter     :: v1_ =1, v2_ =2, v3_ =3
  integer         , parameter     :: nd1_=1, nd2_=2, nd3_=3
  integer         , parameter     :: elm_=1, edg_=2
  double precision, parameter     :: eps         = 1.d-10
  double precision, parameter     :: onethird    = 1.d0 / 3.d0
  double precision, parameter     :: LARGE_VALUE = 1.d20

  allocate( distance(nPoints) )
  distance(:) = LARGE_VALUE

  do iP=1, nPoints

     enclose(:,iP)  = -1
     distance(iP) = LARGE_VALUE

     do iE=1, nElems

        ! ------------------------------------------------------ !
        ! --- [1] edge vector                                --- !
        ! ------------------------------------------------------ !
        edge(x_:z_,e1_) = nodes(x_:z_,elems(nd1_,iE)) - nodes(x_:z_,elems(nd3_,iE) )
        edge(x_:z_,e2_) = nodes(x_:z_,elems(nd2_,iE)) - nodes(x_:z_,elems(nd1_,iE) )
        edge(x_:z_,e3_) = nodes(x_:z_,elems(nd3_,iE)) - nodes(x_:z_,elems(nd2_,iE) )

        ! ------------------------------------------------------ !
        ! --- [2] vector from a vertex to the Point          --- !
        ! ------------------------------------------------------ !
        VtoP(x_:z_,v1_)  = points(x_:z_,iP) - nodes(x_:z_,elems(nd1_,iE) )
        VtoP(x_:z_,v2_)  = points(x_:z_,iP) - nodes(x_:z_,elems(nd2_,iE) )
        VtoP(x_:z_,v3_)  = points(x_:z_,iP) - nodes(x_:z_,elems(nd3_,iE) )

        ! ------------------------------------------------------ !
        ! --- [3] center of gravity / distance from point    --- !
        ! ------------------------------------------------------ !
        cog(x_:z_)       = onethird*( + nodes(x_:z_,elems(nd1_,iE)) &
             &                        + nodes(x_:z_,elems(nd2_,iE)) &
             &                        + nodes(x_:z_,elems(nd3_,iE)) )
        dist             = sqrt( sum( ( points(:,iP) - cog(:) )**2 ) )
        VtoG(x_:z_,v1_)  = cog(x_:z_) - nodes(x_:z_,elems(nd1_,iE) )
        VtoG(x_:z_,v2_)  = cog(x_:z_) - nodes(x_:z_,elems(nd2_,iE) )
        VtoG(x_:z_,v3_)  = cog(x_:z_) - nodes(x_:z_,elems(nd3_,iE) )

        ! ------------------------------------------------------ !
        ! --- [3] cross product                              --- !
        ! ------------------------------------------------------ !
        do iv=v1_, v3_
           crossP(x_,iv) = edge(y_,iv)*VtoP(z_,iv) - edge(z_,iv)*VtoP(y_,iv)
           crossP(y_,iv) = edge(z_,iv)*VtoP(x_,iv) - edge(x_,iv)*VtoP(z_,iv)
           crossP(z_,iv) = edge(x_,iv)*VtoP(y_,iv) - edge(y_,iv)*VtoP(x_,iv)
           crossG(x_,iv) = edge(y_,iv)*VtoG(z_,iv) - edge(z_,iv)*VtoG(y_,iv)
           crossG(y_,iv) = edge(z_,iv)*VtoG(x_,iv) - edge(x_,iv)*VtoG(z_,iv)
           crossG(z_,iv) = edge(x_,iv)*VtoG(y_,iv) - edge(y_,iv)*VtoG(x_,iv)
           normc    (iv) = sqrt( crossP(x_,iv)**2 + crossP(y_,iv)**2 + crossP(z_,iv)**2 )
        enddo

        ! ------------------------------------------------------ !
        ! --- [5] check on edge / on vertex                  --- !
        ! ------------------------------------------------------ !
        do iv=v1_, v3_
           if ( normc(iv) <= eps ) then
              norm_length = ( sqrt( sum( VtoP(:,iv)**2 ) ) ) / sqrt( sum( edge(:,iv)**2 ) )
              if ( ( norm_length >= 0.d0 ).and.( norm_length <= 1.d0 ) ) then
                 ! -- [FOUND] :: exit -- !
                 enclose(elm_,iP) = iE
                 enclose(edg_,iP) = -1
                 exit
              endif
           endif
        enddo

        ! ------------------------------------------------------ !
        ! --- [6] check inside / outside                     --- !
        ! ------------------------------------------------------ !
        do iv=v1_, v3_
           inner(iv) = sum( crossP(:,iv)*crossG(:,iv) )
        enddo
        if ( ( inner(v1_) > 0.d0 ).and.( inner(v2_) > 0.d0 ).and.( inner(v3_) > 0.d0 ) ) then
           ! -- [FOUND] :: exit -- !
           enclose(elm_,iP) = iE
           enclose(edg_,iP) = -1
           exit
        else
           if ( dist < distance(iP) ) then
              do iv=v1_, v3_
                 if ( inner(iv) <= 0.d0  ) then
                    distance(iP)     = min( distance(iP), dist )
                    enclose(elm_,iP) = - iE
                    enclose(edg_,iP) = iv
                 endif
              enddo
           endif
        endif
     enddo

  enddo

  return
end subroutine investigate__inside_triangle


! ====================================================== !
! === interpolation in barycentric coordinate        === !
! ====================================================== !
subroutine barycentric__interpolate( v1, v2, v3, xp )
  implicit none
  double precision, intent(in)    :: v1(3), v2(3), v3(3)
  double precision, intent(inout) :: xp(3)
  double precision                :: denom, pI, w1, w2, w3
  integer         , parameter     :: x_=1, y_=2, z_=3

  denom = ( ( v2(y_)-v3(y_) )*( v1(x_)-v3(x_) ) + ( v3(x_)-v2(x_) )*( v1(y_)-v3(y_) ) )
  if ( denom.eq.0.d0 ) then
     pI = 0.d0
  else
     pI = 1.d0 / denom
  endif
  w1    = pI * ( ( v2(y_)-v3(y_) )*( xp(x_)-v3(x_) ) + ( v3(x_)-v2(x_) )*( xp(y_)-v3(y_) ) )
  w2    = pI * ( ( v3(y_)-v1(y_) )*( xp(x_)-v3(x_) ) + ( v1(x_)-v3(x_) )*( xp(y_)-v3(y_) ) )
  w3    = 1.d0  - w1 - w2

  xp(z_) = w1*v1(z_) + w2*v2(z_) + w3*v3(z_)

  return
end subroutine barycentric__interpolate



! ====================================================== !
! === interpolation from triangle                    === !
! ====================================================== !
subroutine interpolate__trianglePoint( triangle, point )
  implicit none
  double precision, intent(inout) :: triangle(3,3), point(3)
  integer         , parameter     :: x_ =1, y_ =2, z_ =3
  integer         , parameter     :: nSize = 3
  integer                         :: iv
  double precision                :: Amat(3,3), cvec(3), bvec(3)

  do iv=1, nSize
     Amat(iv,1) = 1.d0
     Amat(iv,2) = triangle( x_, iv )
     Amat(iv,3) = triangle( y_, iv )
     bvec(iv)   = triangle( z_, iv )
     cvec(iv)   = 0.d0
  enddo
  call gaussElimin( Amat, cvec, bvec, nSize )
  point(z_)     = cvec(1) + cvec(2)*point(x_) + cvec(3)*point(y_)
  return
end subroutine interpolate__trianglePoint



! ========================================================== !
! === Gauss Elimination Solver                           === !
! ========================================================== !
subroutine gaussElimin( Amat, xvec, bvec, nSize )
  implicit none
  integer         , intent(in)  :: nSize
  double precision, intent(in)  :: Amat(nSize,nSize)
  double precision, intent(in)  :: bvec(nSize)
  double precision, intent(out) :: xvec(nSize)
  integer                       :: i, j, k, ipivot
  double precision              :: Dinv, buff, vpivot
  double precision, parameter   :: eps = 1.d-10
  double precision, allocatable :: Umat(:,:)

  ! ----------------------------------------- !
  ! --- [1] Preparation                   --- !
  ! ----------------------------------------- !
  !  -- [1-1] allocate Umat               --  !
  allocate( Umat(nSize,nSize) )
  !  -- [1-2] Copy AMat & bvec            --  !
  Umat(:,:) = Amat(:,:)
  xvec(:)   = bvec(:)

  ! ----------------------------------------- !
  ! --- [2] Forward Ellimination          --- !
  ! ----------------------------------------- !
  do k=1, nSize

     !  -- [2-1] Pivoting                 --  !
     vpivot = abs( Umat(k,k) )
     ipivot = k
     do j=k+1, nSize
        if ( abs( Umat(j,k) ).gt.vpivot ) then
           vpivot = abs( Umat(j,k) )
           ipivot = j
        endif
     end do
     if ( ipivot.ne.k ) then
        do j=k, nSize
           buff           = Umat(ipivot,j)
           Umat(ipivot,j) = Umat(k     ,j)
           Umat(k     ,j) = buff
        enddo
        buff         = xvec(ipivot)
        xvec(ipivot) = xvec(k)
        xvec(k)      = buff
     end if
     if ( abs( Umat(k,k) ).lt.eps ) then
        write(6,*) '[gaussElimin] Amat :: Singular Matrix :: No Solution End :: @ k= ', k
        stop
     endif
     !  -- [2-2] Diagonal Component       --  !
     Dinv      = 1.d0 / Umat(k,k)
     Umat(k,k) = 1.d0

     !  -- [2-3] Non-Diagonal Component   --  !
     if ( k.eq.nSize ) then
        ! -- [    Last Row :: k == nSize ] -- !
        xvec(k) = Dinv * xvec(k)
     else
        ! -- [Not Last Row :: k != nSize ] -- !
        !  - Division    -  !
        Umat(k,k+1:nSize) = Dinv * Umat(k,k+1:nSize)
        xvec(k)           = Dinv * xvec(k)
        !  - subtraction -  !
        do j=k+1,nSize
           Umat(j,k+1:nSize) = Umat(j,k+1:nSize) - Umat(j,k) * Umat(k,k+1:nSize)
           xvec(j)           = xvec(j)           - Umat(j,k) * xvec(k)
           Umat(j,k)         = 0.d0
        enddo
     endif

  end do

  ! ----------------------------------------- !
  ! --- [3] Backward Substituition        --- !
  ! ----------------------------------------- !
  do k=nSize-1, 1, -1
     do i=nSize, k+1, -1
        xvec(k) = xvec(k) - Umat(k,i)*xvec(i)
     enddo
  enddo

  return
end subroutine gaussElimin
