module difterm_module
  implicit none
  double precision, parameter :: twoThirds = 2.d0/3.d0
  private
  public :: difterm
contains

  subroutine difterm(lo,hi,domlo,domhi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR, NSPEC, QFVAR, QU, QV, QW
    use polyinterp_module, only : cc2zgauss_3d, cc2DzGauss_3d, cc2xyGauss_3d, &
         cc2xyGaussDx_3d, cc2xyGaussDy_3d
    use convert_module, only : cellavg2cc_3d
    use variables_module, only : ctoprim
    use transport_properties, only : get_transport_properties

    integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3), Ulo(3), Uhi(3), &
         fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    integer :: n 
    integer :: g2lo(3), g2hi(3) 
    double precision, allocatable :: Qcc(:,:,:,:)
    double precision, allocatable :: mucc(:,:,:),xicc(:,:,:),lamcc(:,:,:),Ddiacc(:,:,:,:)
    integer :: QzGlo(3), QzGhi(3)
    double precision, allocatable, dimension(:,:,:,:) :: Qz1, Qz2, dveldz1, dveldz2
    double precision, allocatable, dimension(:,:,:) :: muz1, muz2, xiz1, xiz2, lamz1, lamz2
    double precision, allocatable, dimension(:,:,:,:) :: Ddiaz1, Ddiaz2
    integer :: Qzclo(3), Qzchi(3)
    double precision, allocatable, dimension(:,:,:,:) :: Qzc1, Qzc2, Qzc3, Qzc4
    double precision, allocatable, dimension(:,:,:)   :: muzc1,  muzc2,  muzc3,  muzc4
    double precision, allocatable, dimension(:,:,:)   :: xizc1,  xizc2,  xizc3,  xizc4
    double precision, allocatable, dimension(:,:,:)   :: lamzc1, lamzc2, lamzc3, lamzc4
    double precision, allocatable, dimension(:,:,:,:) :: Ddiazc1,Ddiazc2,Ddiazc3,Ddiazc4
    double precision, allocatable, dimension(:,:,:,:) :: dveldx1, dveldx2, dveldx3, dveldx4
    double precision, allocatable, dimension(:,:,:,:) :: dveldy1, dveldy2, dveldy3, dveldy4
    
    g2lo = lo-2 ;    g2hi = hi+2

    allocate(Qcc   (g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3),QFVAR))
    allocate(mucc  (g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3)))
    allocate(xicc  (g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3)))
    allocate(lamcc (g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3)))
    allocate(Ddiacc(g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3),NSPEC))

    do n=1,NVAR
       call cellavg2cc_3d(g2lo,g2hi, U(:,:,:,n), Ulo,Uhi, Qcc(:,:,:,n), g2lo,g2hi)
    end do

    call ctoprim(g2lo, g2hi, Qcc, g2lo, g2hi, QFVAR)

    ! transport coefficients at cell centers
    call get_transport_properties(g2lo,g2hi, Qcc,g2lo,g2hi,QFVAR, &
         mucc,xicc,lamcc,Ddiacc, g2lo,g2hi)

    QzGlo(1:2) = g2lo(1:2)
    QzGlo(3) = lo(3)
    QzGhi(1:2) = g2hi(1:2)
    QzGhi(3) = hi(3)

    allocate(   Qz1(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3),QFVAR))
    allocate(  muz1(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3)))
    allocate(  xiz1(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3)))
    allocate( lamz1(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3)))
    allocate(Ddiaz1(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3),NSPEC))

    allocate(   Qz2(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3),QFVAR))
    allocate(  muz2(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3)))
    allocate(  xiz2(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3)))
    allocate( lamz2(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3)))
    allocate(Ddiaz2(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3),NSPEC))

    allocate(dveldz1(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3),3))
    allocate(dveldz2(QzGlo(1):QzGhi(1),QzGlo(2):QzGhi(2),QzGlo(3):QzGhi(3),3))

    call cc2zgauss_3d(QzGlo,QzGhi,  mucc, g2lo,g2hi,  muz1,  muz2, QzGlo,QzGhi)
    call cc2zgauss_3d(QzGlo,QzGhi,  xicc, g2lo,g2hi,  xiz1,  xiz2, QzGlo,QzGhi)
    call cc2zgauss_3d(QzGlo,QzGhi, lamcc, g2lo,g2hi, lamz1, lamz2, QzGlo,QzGhi)
    do n=1,NSPEC
       call cc2zgauss_3d(QzGlo,QzGhi, Ddiacc(:,:,:,n), g2lo,g2hi, &
            Ddiaz1(:,:,:,n), Ddiaz2(:,:,:,n), QzGlo,QzGhi)
    end do

    do n=1,QFVAR
       call cc2zgauss_3d(QzGlo,QzGhi, Qcc(:,:,:,n), g2lo,g2hi, &
            Qz1(:,:,:,n), Qz2(:,:,:,n), QzGlo,QzGhi)
    end do

    do n=1,3
       call cc2DzGauss_3d(QzGlo,QzGhi, Qcc(:,:,:,QU+n-1), g2lo,g2hi, &
            dveldz1(:,:,:,n), dveldz2(:,:,:,n), QzGlo,QzGhi)
    end do

    call diff_xy(lo,hi, domlo, domhi, dveldz1, Qz1, muz1, xiz1, lamz1, Ddiaz1, QzGlo, QzGhi, &
         fx, fxlo, fxhi, fy, fylo, fyhi, dx)

    call diff_xy(lo,hi, domlo, domhi, dveldz2, Qz2, muz2, xiz2, lamz2, Ddiaz2, QzGlo, QzGhi, &
         fx, fxlo, fxhi, fy, fylo, fyhi, dx)

    deallocate(Qz1,Qz2,dveldz1,dveldz2)
    deallocate(muz1,muz2,xiz1,xiz2,lamz1,lamz2,Ddiaz1,Ddiaz2)

    ! work on z-direction

    Qzclo(1:2) = lo(1:2)
    Qzclo(3) = lo(3)-2
    Qzchi(1:2) = hi(1:2)
    Qzchi(3) = hi(3)+2

    allocate(Qzc1(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),QFVAR))
    allocate(Qzc2(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),QFVAR))
    allocate(Qzc3(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),QFVAR))
    allocate(Qzc4(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),QFVAR))

    allocate(muzc1(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))
    allocate(muzc2(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))
    allocate(muzc3(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))
    allocate(muzc4(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))

    allocate(xizc1(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))
    allocate(xizc2(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))
    allocate(xizc3(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))
    allocate(xizc4(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))

    allocate(lamzc1(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))
    allocate(lamzc2(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))
    allocate(lamzc3(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))
    allocate(lamzc4(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3)))

    allocate(Ddiazc1(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),nspec))
    allocate(Ddiazc2(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),nspec))
    allocate(Ddiazc3(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),nspec))
    allocate(Ddiazc4(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),nspec))

    allocate(dveldx1(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),2))
    allocate(dveldx2(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),2))
    allocate(dveldx3(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),2))
    allocate(dveldx4(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),2))

    allocate(dveldy1(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),2))
    allocate(dveldy2(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),2))
    allocate(dveldy3(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),2))
    allocate(dveldy4(Qzclo(1):Qzchi(1),Qzclo(2):Qzchi(2),Qzclo(3):Qzchi(3),2))

    do n=1,QFVAR
       call cc2xyGauss_3d(Qzclo,Qzchi,Qcc(:,:,:,n),g2lo,g2hi, &
            Qzc1(:,:,:,n),Qzc2(:,:,:,n),Qzc3(:,:,:,n),Qzc4(:,:,:,n),Qzclo,Qzchi)
    end do
    call cc2xyGauss_3d(Qzclo,Qzchi, mucc,g2lo,g2hi, muzc1, muzc2, muzc3, muzc4,Qzclo,Qzchi)
    call cc2xyGauss_3d(Qzclo,Qzchi, xicc,g2lo,g2hi, xizc1, xizc2, xizc3, xizc4,Qzclo,Qzchi)
    call cc2xyGauss_3d(Qzclo,Qzchi,lamcc,g2lo,g2hi,lamzc1,lamzc2,lamzc3,lamzc4,Qzclo,Qzchi)
    do n=1,nspec
       call cc2xyGauss_3d(Qzclo,Qzchi,Ddiacc(:,:,:,n),g2lo,g2hi, &
            Ddiazc1(:,:,:,n),Ddiazc2(:,:,:,n),Ddiazc3(:,:,:,n),Ddiazc4(:,:,:,n),Qzclo,Qzchi)
    end do

    ! du/dx
    call cc2xyGaussDx_3d(Qzclo,Qzchi,Qcc(:,:,:,QU),g2lo,g2hi, &
         dveldx1(:,:,:,1),dveldx2(:,:,:,1),dveldx3(:,:,:,1),dveldx4(:,:,:,1),Qzclo,Qzchi)
    ! dw/dx
    call cc2xyGaussDx_3d(Qzclo,Qzchi,Qcc(:,:,:,QW),g2lo,g2hi, &
         dveldx1(:,:,:,2),dveldx2(:,:,:,2),dveldx3(:,:,:,2),dveldx4(:,:,:,2),Qzclo,Qzchi)
    
    ! dv/dy
    call cc2xyGaussDy_3d(Qzclo,Qzchi,Qcc(:,:,:,QV),g2lo,g2hi, &
         dveldy1(:,:,:,1),dveldy2(:,:,:,1),dveldy3(:,:,:,1),dveldy4(:,:,:,1),Qzclo,Qzchi)
    ! dw/dy
    call cc2xyGaussDy_3d(Qzclo,Qzchi,Qcc(:,:,:,QW),g2lo,g2hi, &
         dveldy1(:,:,:,2),dveldy2(:,:,:,2),dveldy3(:,:,:,2),dveldy4(:,:,:,2),Qzclo,Qzchi)

    deallocate(Qcc,mucc,xicc,lamcc,Ddiacc)

    call diff_z(lo,hi,domlo,domhi,Qzc1,muzc1,xizc1,lamzc1,Ddiazc1,dveldx1,dveldy1,Qzclo,Qzchi, &
         fz, fzlo, fzhi, dx)
    call diff_z(lo,hi,domlo,domhi,Qzc2,muzc2,xizc2,lamzc2,Ddiazc2,dveldx2,dveldy2,Qzclo,Qzchi, &
         fz, fzlo, fzhi, dx)
    call diff_z(lo,hi,domlo,domhi,Qzc3,muzc3,xizc3,lamzc3,Ddiazc3,dveldx3,dveldy3,Qzclo,Qzchi, &
         fz, fzlo, fzhi, dx)
    call diff_z(lo,hi,domlo,domhi,Qzc4,muzc4,xizc4,lamzc4,Ddiazc4,dveldx4,dveldy4,Qzclo,Qzchi, &
         fz, fzlo, fzhi, dx)

    deallocate(Qzc1,Qzc2,Qzc3,Qzc4)
    deallocate(muzc1,  muzc2,  muzc3,  muzc4)     
    deallocate(xizc1,  xizc2,  xizc3,  xizc4)     
    deallocate(lamzc1, lamzc2, lamzc3, lamzc4)  
    deallocate(Ddiazc1,Ddiazc2,Ddiazc3,Ddiazc4)
    deallocate(dveldx1, dveldx2, dveldx3, dveldx4)
    deallocate(dveldy1, dveldy2, dveldy3, dveldy4)

    return
  end subroutine difterm


  subroutine diff_xy(lo, hi, domlo, domhi, dveldz, Q, mu, xi, lam, Ddia, qlo, qhi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, dx)
    
    use meth_params_module, only : NVAR, NSPEC, QCVAR, QFVAR, QU, QV
    use polyinterp_module, only : cc2xface_2d, cc2yface_2d, cc2DxYface_2d, cc2DyXface_2d

    integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3), qlo(3), qhi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::dveldz( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),3)
    double precision,intent(in   )::     Q( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QFVAR) 
    double precision,intent(in   )::    mu( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3)) 
    double precision,intent(in   )::    xi( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3)) 
    double precision,intent(in   )::   lam( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3)) 
    double precision,intent(in   )::  Ddia( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),NSPEC) 
    double precision,intent(inout)::    fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::    fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)

    double precision :: dxinv(3)
    integer :: k, n
    integer :: g2lo(2), g2hi(2)
    double precision, allocatable :: Qc1(:,:,:), Qc2(:,:,:)
    double precision, allocatable :: tmp1(:,:), tmp2(:,:)
    integer :: Qflo(2), Qfhi(2)
    double precision, allocatable :: Qf1(:,:,:), Qf2(:,:,:)
    double precision, dimension(:,:), allocatable :: mu1, mu2, xi1, xi2, lam1, lam2
    double precision, dimension(:,:,:), allocatable :: Ddia1, Ddia2
    double precision, allocatable :: dveldz1(:,:,:), dveldz2(:,:,:)
    double precision, allocatable :: dveldt1(:,:,:), dveldt2(:,:,:)  ! dt means dx or dy
    integer :: tlo(2), thi(2)
    double precision, parameter :: fac = 0.25d0 ! due to Gauss quadrature

    dxinv = 1.d0/dx

    g2lo = lo(1:2)-2  
    g2hi = hi(1:2)+2

    allocate(Qc1(g2lo(1):g2hi(1),g2lo(2):g2hi(2),QCVAR))
    allocate(Qc2(g2lo(1):g2hi(1),g2lo(2):g2hi(2),QCVAR))    

    allocate(tmp1(g2lo(1):g2hi(1),g2lo(2):g2hi(2)))
    allocate(tmp2(g2lo(1):g2hi(1),g2lo(2):g2hi(2)))    

    Qflo = lo(1:2)
    Qfhi = hi(1:2)+1

    allocate(dveldt1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2))
    allocate(dveldt2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2))

    allocate(Qf1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR))
    allocate(Qf2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR))

    allocate(dveldz1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2))
    allocate(dveldz2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2))

    allocate(  mu1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(  xi1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate( lam1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(Ddia1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC))

    allocate(  mu2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(  xi2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate( lam2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(Ddia2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC))

    do k = lo(3), hi(3)

       ! ----- compute x-direction flux first -----

       ! cell center => Gauss points on x-face
       call cc2xface_2d(lo(1:2),hi(1:2), mu(:,:,k),qlo(1:2),qhi(1:2), mu1, mu2,Qflo,Qfhi,&
            tmp1, tmp2, g2lo, g2hi)
       call cc2xface_2d(lo(1:2),hi(1:2), xi(:,:,k),qlo(1:2),qhi(1:2), xi1, xi2,Qflo,Qfhi,&
            tmp1, tmp2, g2lo, g2hi)
       call cc2xface_2d(lo(1:2),hi(1:2),lam(:,:,k),qlo(1:2),qhi(1:2),lam1,lam2,Qflo,Qfhi,&
            tmp1, tmp2, g2lo, g2hi)
       do n=1,NSPEC
          call cc2xface_2d(lo(1:2),hi(1:2),Ddia(:,:,k,n),qlo(1:2),qhi(1:2), &
               Ddia1(:,:,n),Ddia2(:,:,n),Qflo,Qfhi, tmp1, tmp2, g2lo, g2hi)
       end do

       ! dU/dz
       call cc2xface_2d(lo(1:2),hi(1:2),dveldz(:,:,k,1),qlo(1:2),qhi(1:2), &
            dveldz1(:,:,1),dveldz2(:,:,1),Qflo,Qfhi, tmp1, tmp2, g2lo, g2hi)
       ! dW/dz
       call cc2xface_2d(lo(1:2),hi(1:2),dveldz(:,:,k,3),qlo(1:2),qhi(1:2), &
            dveldz1(:,:,2),dveldz2(:,:,2),Qflo,Qfhi, tmp1, tmp2, g2lo, g2hi)

       ! cell-center => Qc: center-in-x and Gauss-point-in-y 
       !                Qf: xface and Gauss-point-in-y
       do n=1,QCVAR
          call cc2xface_2d(lo(1:2),hi(1:2), Q(:,:,k,n), qlo(1:2),qhi(1:2), &
               Qf1(:,:,n), Qf2(:,:,n), Qflo, Qfhi, &
               Qc1(:,:,n), Qc2(:,:,n), g2lo, g2hi)
       end do
       do n=QCVAR+1,QFVAR
          call cc2xface_2d(lo(1:2),hi(1:2), Q(:,:,k,n), qlo(1:2),qhi(1:2), &
               Qf1(:,:,n), Qf2(:,:,n), Qflo, Qfhi, tmp1, tmp2, g2lo, g2hi)
       end do

       ! cell-average of ? => xface & Gauss-point-in-y of d?/dy
       call cc2DyXface_2d(lo(1:2),hi(1:2), Q(:,:,k,QU), qlo(1:2), qhi(1:2), &
            dveldt1(:,:,1), dveldt2(:,:,1), Qflo, Qfhi, tmp1, tmp2, g2lo, g2hi)
       call cc2DyXface_2d(lo(1:2),hi(1:2), Q(:,:,k,QV), qlo(1:2), qhi(1:2), &
            dveldt1(:,:,2), dveldt2(:,:,2), Qflo, Qfhi, tmp1, tmp2, g2lo, g2hi)

       tlo = lo(1:2)
       thi(1) = hi(1)+1
       thi(2) = hi(2)
       call comp_diff_flux_x(tlo, thi, k, fx, fxlo, fxhi, &
            Qf1, mu1, xi1, lam1, Ddia1, dveldt1, dveldz1, Qflo, Qfhi, &
            Qc1, g2lo, g2hi, dxinv, fac, domlo, domhi)
       call comp_diff_flux_x(tlo, thi, k, fx, fxlo, fxhi, &
            Qf2, mu2, xi2, lam2, Ddia2, dveldt2, dveldz2, Qflo, Qfhi, &
            Qc2, g2lo, g2hi, dxinv, fac, domlo, domhi)

       ! ----- compute y-direction flux -----

       ! cell center => Gauss points on y-face
       call cc2yface_2d(lo(1:2),hi(1:2), mu(:,:,k),qlo(1:2),qhi(1:2), mu1, mu2,Qflo,Qfhi,&
            tmp1, tmp2, g2lo, g2hi)
       call cc2yface_2d(lo(1:2),hi(1:2), xi(:,:,k),qlo(1:2),qhi(1:2), xi1, xi2,Qflo,Qfhi,&
            tmp1, tmp2, g2lo, g2hi)
       call cc2yface_2d(lo(1:2),hi(1:2),lam(:,:,k),qlo(1:2),qhi(1:2),lam1,lam2,Qflo,Qfhi,&
            tmp1, tmp2, g2lo, g2hi)
       do n=1,NSPEC
          call cc2yface_2d(lo(1:2),hi(1:2),Ddia(:,:,k,n),qlo(1:2),qhi(1:2), &
               Ddia1(:,:,n),Ddia2(:,:,n),Qflo,Qfhi, tmp1, tmp2, g2lo, g2hi)
       end do

       ! dV/dz
       call cc2yface_2d(lo(1:2),hi(1:2),dveldz(:,:,k,2),qlo(1:2),qhi(1:2), &
            dveldz1(:,:,1),dveldz2(:,:,1),Qflo,Qfhi, tmp1, tmp2, g2lo, g2hi)
       ! dW/dz
       call cc2yface_2d(lo(1:2),hi(1:2),dveldz(:,:,k,3),qlo(1:2),qhi(1:2), &
            dveldz1(:,:,2),dveldz2(:,:,2),Qflo,Qfhi, tmp1, tmp2, g2lo, g2hi)

       ! cell-center => Qc: center-in-y and Gauss-point-in-x
       !                Qf: yface and Gauss-point-in-x
       do n=1,QCVAR
          call cc2yface_2d(lo(1:2),hi(1:2), Q(:,:,k,n), qlo(1:2),qhi(1:2), &
               Qf1(:,:,n), Qf2(:,:,n), Qflo, Qfhi, &
               Qc1(:,:,n), Qc2(:,:,n), g2lo, g2hi)
       end do
       do n=QCVAR+1,QFVAR
          call cc2yface_2d(lo(1:2),hi(1:2), Q(:,:,k,n), qlo(1:2),qhi(1:2), &
               Qf1(:,:,n), Qf2(:,:,n), Qflo, Qfhi, tmp1, tmp2, g2lo, g2hi)
       end do

       ! cell-average of ? => yface & Gauss-point-in-x of d?/dx
       call cc2DxYface_2d(lo(1:2),hi(1:2), Q(:,:,k,QU), qlo(1:2), qhi(1:2), &
            dveldt1(:,:,1), dveldt2(:,:,1), Qflo, Qfhi, tmp1, tmp2, g2lo, g2hi)
       call cc2DxYface_2d(lo(1:2),hi(1:2), Q(:,:,k,QV), qlo(1:2), qhi(1:2), &
            dveldt1(:,:,2), dveldt2(:,:,2), Qflo, Qfhi, tmp1, tmp2, g2lo, g2hi)

       tlo = lo(1:2)
       thi(1) = hi(1)
       thi(2) = hi(2)+1
       call comp_diff_flux_y(tlo, thi, k, fy, fylo, fyhi, &
            Qf1, mu1, xi1, lam1, Ddia1, dveldt1, dveldz1, Qflo, Qfhi, &
            Qc1, g2lo, g2hi, dxinv, fac, domlo, domhi)
       call comp_diff_flux_y(tlo, thi, k, fy, fylo, fyhi, &
            Qf2, mu2, xi2, lam2, Ddia2, dveldt2, dveldz2, Qflo, Qfhi, &
            Qc2, g2lo, g2hi, dxinv, fac, domlo, domhi)

    end do

    deallocate(Qc1,Qc2,tmp1,tmp2,dveldt1,dveldt2,Qf1,Qf2,dveldz1,dveldz2)
    deallocate(mu1,xi1,lam1,Ddia1)
    deallocate(mu2,xi2,lam2,Ddia2)

  end subroutine diff_xy


  subroutine comp_diff_flux_x(lo, hi, k, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dveldy, dveldz, Qflo, Qfhi, &
       Qc, Qclo, Qchi, dxinv, fac, domlo, domhi)

    use meth_params_module
    use derivative_stencil_module, only : FD4

    integer, intent(in) :: lo(2), hi(2), k, flo(3), fhi(3), Qflo(2), Qfhi(2), Qclo(2), Qchi(2), &
         domlo(3), domhi(3)
    double precision, intent(in) :: dxinv(3), fac
    double precision, intent(inout) ::   flx( flo(1): fhi(1), flo(2): fhi(2),flo(3):fhi(3),NVAR)
    double precision, intent(in   ) ::    Qf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR)
    double precision, intent(in   ) ::    mu(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::    xi(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::   lam(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::  Ddia(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC)
    double precision, intent(in   ) ::dveldy(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2)
    double precision, intent(in   ) ::dveldz(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2)
    double precision, intent(in   ) ::    Qc(Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR)

    integer :: i, j, n, UYN, QYN, QXN, QHN
    double precision :: tauxx, tauxy, tauxz
    double precision :: dudx, dudy, dudz, dvdx, dvdy, dwdx, dwdz, divu
    double precision :: dTdx, dXdx, Vd
    double precision :: ek, rhovn
    double precision, dimension(lo(1):hi(1)) :: dlnpdx, Vc

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ! viscous stress
          dudx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QU) + FD4(-1)*Qc(i-1,j,QU) &
               + FD4(0)*Qc(i,j,QU) + FD4(1)*Qc(i+1,j,QU))
          dvdx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QV) + FD4(-1)*Qc(i-1,j,QV) &
               + FD4(0)*Qc(i,j,QV) + FD4(1)*Qc(i+1,j,QV))
          dwdx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QW) + FD4(-1)*Qc(i-1,j,QW) &
               + FD4(0)*Qc(i,j,QW) + FD4(1)*Qc(i+1,j,QW))
          dudy = dxinv(2)*dveldy(i,j,1)
          dvdy = dxinv(2)*dveldy(i,j,2)
          dudz = dxinv(3)*dveldz(i,j,1)
          dwdz = dxinv(3)*dveldz(i,j,2)
          divu = dudx + dvdy + dwdz
          tauxx = mu(i,j)*(2.d0*dudx-twoThirds*divu) + xi(i,j)*divu
          tauxy = mu(i,j)*(dudy+dvdx)
          tauxz = mu(i,j)*(dudz+dwdx)
          flx(i,j,k,UMX)   = flx(i,j,k,UMX)   - fac*tauxx
          flx(i,j,k,UMY)   = flx(i,j,k,UMY)   - fac*tauxy
          flx(i,j,k,UMZ)   = flx(i,j,k,UMZ)   - fac*tauxz
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac* &
               (tauxx*Qf(i,j,QU) + tauxy*Qf(i,j,QV) + tauxz*Qf(i,j,QW))

          ! thermal conduction
          dTdx = dxinv(1) * (FD4(-2)*Qc(i-2,j,QTEMP) + FD4(-1)*Qc(i-1,j,QTEMP) &
               + FD4(0)*Qc(i,j,QTEMP) + FD4(1)*Qc(i+1,j,QTEMP))
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac*lam(i,j)*dTdx

          ! compute dpdx
          dlnpdx(i) = dxinv(1) * (FD4(-2)*Qc(i-2,j,QPRES) + FD4(-1)*Qc(i-1,j,QPRES) &
               + FD4(0)*Qc(i,j,QPRES) + FD4(1)*Qc(i+1,j,QPRES)) / Qf(i,j,QPRES)
          Vc(i) = 0.d0
       end do

       do n=1,NSPEC
          UYN = UFS+n-1
          QYN = QFY+n-1
          QXN = QFX+n-1
          QHN = QFH+n-1
          do i = lo(1), hi(1)
             dXdx = dxinv(1) * (FD4(-2)*Qc(i-2,j,QXN) + FD4(-1)*Qc(i-1,j,QXN) &
                  + FD4(0)*Qc(i,j,QXN) + FD4(1)*Qc(i+1,j,QXN))
             Vd = -Ddia(i,j,n)*(dXdx + (Qf(i,j,QXN)-Qf(i,j,QYN))*dlnpdx(i))
             
             flx(i,j,k,UYN) = flx(i,j,k,UYN) + fac*Vd
             Vc(i) = Vc(i) + Vd
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) + fac*Vd*Qf(i,j,QHN)
          end do
       end do

       do n=1,NSPEC
          UYN = UFS+n-1
          QYN = QFY+n-1
          QHN = QFH+n-1
          do i = lo(1), hi(1)
             flx(i,j,k,UYN )  = flx(i,j,k,UYN  ) - (fac*Qf(i,j,QYN)*Vc(i))
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - (fac*Qf(i,j,QYN)*Vc(i))*Qf(i,j,QHN)
          end do
       end do
    end do

    if (.not. do_weno) then
       ! compute hyperbolic flux
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhovn = fac*Qf(i,j,QRHO)*Qf(i,j,QU)
             flx(i,j,k,URHO) = flx(i,j,k,URHO) + rhovn
             flx(i,j,k,UMX ) = flx(i,j,k,UMX ) + rhovn*Qf(i,j,QU) + fac*Qf(i,j,QPRES)
             flx(i,j,k,UMY ) = flx(i,j,k,UMY ) + rhovn*Qf(i,j,QV)
             flx(i,j,k,UMZ ) = flx(i,j,k,UMZ ) + rhovn*Qf(i,j,QW)

             ek = 0.5d0*(Qf(i,j,QU)**2+Qf(i,j,QV)**2+Qf(i,j,QW)**2)
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) + rhovn*ek
             do n=1,NSPEC
                flx(i,j,k,UEDEN)   = flx(i,j,k,UEDEN)   + (rhovn*Qf(i,j,QFY+n-1))*Qf(i,j,QFH+n-1)
                flx(i,j,k,UFS+n-1) = flx(i,j,k,UFS+n-1) + (rhovn*Qf(i,j,QFY+n-1))
             end do
          end do
       end do
    end if

  end subroutine comp_diff_flux_x


  subroutine comp_diff_flux_y(lo, hi, k, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dveldx, dveldz, Qflo, Qfhi, &
       Qc, Qclo, Qchi, dxinv, fac, domlo, domhi)

    use prob_params_module, only : physbc_lo, physbc_hi, NoSlipWall
    use meth_params_module
    use derivative_stencil_module, only : FD4
    use RNS_boundary_module, only : Twall

    integer, intent(in) :: lo(2), hi(2), k, flo(3), fhi(3), Qflo(2), Qfhi(2), Qclo(2), Qchi(2), &
         domlo(3), domhi(3)
    double precision, intent(in) :: dxinv(3), fac
    double precision, intent(inout) ::   flx( flo(1): fhi(1), flo(2): fhi(2),flo(3):fhi(3),NVAR)
    double precision, intent(in   ) ::    Qf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR)
    double precision, intent(in   ) ::    mu(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::    xi(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::   lam(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::  Ddia(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC)
    double precision, intent(in   ) ::dveldx(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2)
    double precision, intent(in   ) ::dveldz(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),2)
    double precision, intent(in   ) ::    Qc(Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR)

    integer :: i, j, n, UYN, QYN, QXN, QHN, jstart
    double precision :: tauyy, tauxy, tauyz
    double precision :: dudx, dudy, dvdx, dvdy, dvdz, dwdy, dwdz, divu
    double precision :: dTdy, dXdy, Vd
    double precision :: ek, rhovn
    double precision :: dlnpdytmp, Vctmp, foo
    double precision, allocatable :: dlnpdy(:,:), Vc(:,:)

    jstart = lo(2)

    if (lo(2).eq.domlo(2) .and. physbc_lo(2).eq.NoSlipWall) then
       jstart = lo(2)+2
       
       ! the face on the wall
       j = lo(2)
       do i=lo(1),hi(1)
          ! viscous stress
          dudy = 0.5d0*dxinv(2)*Qc(i,lo(2),QU)
          dvdy = 0.5d0*dxinv(2)*Qc(i,lo(2),QV)
          dwdy = 0.5d0*dxinv(2)*Qc(i,lo(2),QW)
          dudx = 0.d0
          dvdx = 0.d0
          dvdz = 0.d0
          dwdz = 0.d0
          divu = dudx + dvdy + dwdz
          tauyy = mu(i,j)*(2.d0*dvdy-twoThirds*divu) + xi(i,j)*divu
          tauxy = mu(i,j)*(dudy+dvdx)
          tauyz = mu(i,j)*(dwdy+dvdz)
          flx(i,j,k,UMX)   = flx(i,j,k,UMX)   - fac*tauxy
          flx(i,j,k,UMY)   = flx(i,j,k,UMY)   - fac*tauyy
          flx(i,j,k,UMZ)   = flx(i,j,k,UMZ)   - fac*tauyz

          ! thermal conduction
          dTdy = 0.5d0*dxinv(2)*(Qc(i,lo(2),QTEMP)-Twall)
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac*lam(i,j)*dTdy
       end do

       ! the next face
       j = lo(2)+1
       do i=lo(1),hi(1)
          ! viscous stress
          dudy = dxinv(2)*(Qc(i,lo(2)+1,QU)-Qc(i,lo(2),QU))
          dvdy = dxinv(2)*(Qc(i,lo(2)+1,QV)-Qc(i,lo(2),QV))
          dwdy = dxinv(2)*(Qc(i,lo(2)+1,QW)-Qc(i,lo(2),QW))
          dudx = dxinv(1)*dveldx(i,j,1)
          dvdx = dxinv(1)*dveldx(i,j,2)
          dvdz = dxinv(3)*dveldz(i,j,1)
          dwdz = dxinv(3)*dveldz(i,j,2)
          divu = dudx + dvdy + dwdz
          tauyy = mu(i,j)*(2.d0*dvdy-twoThirds*divu) + xi(i,j)*divu
          tauxy = mu(i,j)*(dudy+dvdx)
          tauyz = mu(i,j)*(dwdy+dvdz)
          flx(i,j,k,UMX)   = flx(i,j,k,UMX)   - fac*tauxy
          flx(i,j,k,UMY)   = flx(i,j,k,UMY)   - fac*tauyy
          flx(i,j,k,UMZ)   = flx(i,j,k,UMZ)   - fac*tauyz
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac* &
               0.5d0*( tauxy*(Qc(i,lo(2),QU)+Qc(i,lo(2)+1,QU)) &
               &     + tauyy*(Qc(i,lo(2),QV)+Qc(i,lo(2)+1,QV)) &
               &     + tauyz*(Qc(i,lo(2),QW)+Qc(i,lo(2)+1,QW)) )

          ! thermal conduction
          dTdy = dxinv(2) * (Qc(i,lo(2)+1,QTEMP) - Qc(i,lo(2),QTEMP))
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac*lam(i,j)*dTdy

          ! compute dpdy
          dlnpdytmp = dxinv(2) * (Qc(i,lo(2)+1,QPRES) - Qc(i,lo(2),QPRES)) &
               / (0.5d0*(Qc(i,lo(2),QPRES) + Qc(i,lo(2)+1,QPRES)))
          Vctmp = 0.d0

          do n=1,NSPEC
             UYN = UFS+n-1
             QYN = QFY+n-1
             QXN = QFX+n-1
             QHN = QFH+n-1

             dXdy = dxinv(2) * (Qc(i,lo(2)+1,QXN) - Qc(i,lo(2),QXN))
             Vd = -Ddia(i,j,n)*(dXdy + dlnpdytmp* &
                  0.5d0*(Qc(i,lo(2),QXN)+Qc(i,lo(2)+1,QXN) &
                  &    - Qc(i,lo(2),QYN)-Qc(i,lo(2)+1,QYN)) ) 
             flx(i,j,k,UYN) = flx(i,j,k,UYN) + fac*Vd
             Vctmp = Vctmp + Vd
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) + fac*Vd*Qf(i,j,QHN)
          end do

          do n=1,NSPEC
             UYN = UFS+n-1
             QYN = QFY+n-1
             QHN = QFH+n-1
             foo = fac*0.5d0*(Qc(i,lo(2),QYN)+Qc(i,lo(2)+1,QYN))*Vctmp
             flx(i,j,k,UYN )  = flx(i,j,k,UYN  ) - foo
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - foo*Qf(i,j,QHN)
          end do
       end do
    end if

    allocate(dlnpdy(lo(1):hi(1),jstart:hi(2)))
    allocate(    Vc(lo(1):hi(1),jstart:hi(2)))

    do j=jstart,hi(2)
       do i=lo(1),hi(1)
          ! viscous stress
          dudy = dxinv(2)*(FD4(-2)*Qc(i,j-2,QU) + FD4(-1)*Qc(i,j-1,QU) &
               + FD4(0)*Qc(i,j,QU) + FD4(1)*Qc(i,j+1,QU))
          dvdy = dxinv(2)*(FD4(-2)*Qc(i,j-2,QV) + FD4(-1)*Qc(i,j-1,QV) &
               + FD4(0)*Qc(i,j,QV) + FD4(1)*Qc(i,j+1,QV))
          dwdy = dxinv(2)*(FD4(-2)*Qc(i,j-2,QW) + FD4(-1)*Qc(i,j-1,QW) &
               + FD4(0)*Qc(i,j,QW) + FD4(1)*Qc(i,j+1,QW))
          dudx = dxinv(1)*dveldx(i,j,1)
          dvdx = dxinv(1)*dveldx(i,j,2)
          dvdz = dxinv(3)*dveldz(i,j,1)
          dwdz = dxinv(3)*dveldz(i,j,2)
          divu = dudx + dvdy + dwdz
          tauyy = mu(i,j)*(2.d0*dvdy-twoThirds*divu) + xi(i,j)*divu
          tauxy = mu(i,j)*(dudy+dvdx)
          tauyz = mu(i,j)*(dwdy+dvdz)
          flx(i,j,k,UMX)   = flx(i,j,k,UMX)   - fac*tauxy
          flx(i,j,k,UMY)   = flx(i,j,k,UMY)   - fac*tauyy
          flx(i,j,k,UMZ)   = flx(i,j,k,UMZ)   - fac*tauyz
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac* &
               (tauxy*Qf(i,j,QU) + tauyy*Qf(i,j,QV) + tauyz*Qf(i,j,QW))

          ! thermal conduction
          dTdy = dxinv(2) * (FD4(-2)*Qc(i,j-2,QTEMP) + FD4(-1)*Qc(i,j-1,QTEMP) &
               + FD4(0)*Qc(i,j,QTEMP) + FD4(1)*Qc(i,j+1,QTEMP))
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac*lam(i,j)*dTdy

          ! compute dpdy
          dlnpdy(i,j) = dxinv(2) * (FD4(-2)*Qc(i,j-2,QPRES) + FD4(-1)*Qc(i,j-1,QPRES) &
               + FD4(0)*Qc(i,j,QPRES) + FD4(1)*Qc(i,j+1,QPRES)) / Qf(i,j,QPRES)
          Vc(i,j) = 0.d0
       end do
    end do

    do n=1,NSPEC
       UYN = UFS+n-1
       QYN = QFY+n-1
       QXN = QFX+n-1
       QHN = QFH+n-1
       do j = jstart, hi(2)
          do i = lo(1), hi(1)
             dXdy = dxinv(2) * (FD4(-2)*Qc(i,j-2,QXN) + FD4(-1)*Qc(i,j-1,QXN) &
                  + FD4(0)*Qc(i,j,QXN) + FD4(1)*Qc(i,j+1,QXN))
             Vd = -Ddia(i,j,n)*(dXdy + (Qf(i,j,QXN)-Qf(i,j,QYN))*dlnpdy(i,j))
             
             flx(i,j,k,UYN) = flx(i,j,k,UYN) + fac*Vd
             Vc(i,j) = Vc(i,j) + Vd
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) + fac*Vd*Qf(i,j,QHN)
          end do
       end do
    end do

    do n=1,NSPEC
       UYN = UFS+n-1
       QYN = QFY+n-1
       QHN = QFH+n-1
       do j= jstart, hi(2)
          do i = lo(1), hi(1)
             flx(i,j,k,UYN )  = flx(i,j,k,UYN  ) - (fac*Qf(i,j,QYN)*Vc(i,j))
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - (fac*Qf(i,j,QYN)*Vc(i,j))*Qf(i,j,QHN)
          end do
       end do
    end do

    deallocate(dlnpdy,Vc)

    if (.not. do_weno) then
       ! compute hyperbolic flux
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhovn = fac*Qf(i,j,QRHO)*Qf(i,j,QV)
             flx(i,j,k,URHO) = flx(i,j,k,URHO) + rhovn
             flx(i,j,k,UMX ) = flx(i,j,k,UMX ) + rhovn*Qf(i,j,QU)
             flx(i,j,k,UMY ) = flx(i,j,k,UMY ) + rhovn*Qf(i,j,QV) + fac*Qf(i,j,QPRES)
             flx(i,j,k,UMZ ) = flx(i,j,k,UMZ ) + rhovn*Qf(i,j,QW)

             ek = 0.5d0*(Qf(i,j,QU)**2+Qf(i,j,QV)**2+Qf(i,j,QW)**2)
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) + rhovn*ek
             do n=1,NSPEC
                flx(i,j,k,UEDEN)   = flx(i,j,k,UEDEN)   + (rhovn*Qf(i,j,QFY+n-1))*Qf(i,j,QFH+n-1)
                flx(i,j,k,UFS+n-1) = flx(i,j,k,UFS+n-1) + (rhovn*Qf(i,j,QFY+n-1))
             end do
          end do
       end do
    end if

  end subroutine comp_diff_flux_y


  subroutine diff_z(lo, hi, domlo, domhi, Q, mu, xi, lam, Ddia, dveldx, dveldy, qlo, qhi, &
       f, flo, fhi, dx)
    
    use meth_params_module, only : NVAR, NSPEC, QFVAR
    use polyinterp_module, only : cc2zface_3d
    
    integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3), qlo(3), qhi(3), flo(3), fhi(3)
    double precision, intent(in   )::     Q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),QFVAR)
    double precision, intent(in   )::    mu(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),3)
    double precision, intent(in   )::    xi(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(in   )::   lam(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(in   )::  Ddia(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),NSPEC)
    double precision, intent(in   )::dveldx(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),2)
    double precision, intent(in   )::dveldy(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),2)
    double precision, intent(inout)::     f(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),NVAR)
    double precision, intent(in) :: dx(3)

    double precision  :: dxinv(3)
    integer :: n, Qflo(3), Qfhi(3)
    double precision, dimension(:,:,:,:), allocatable :: Qf, dvdxf, dvdyf, Ddiaf
    double precision, dimension(:,:,:)  , allocatable :: muf, xif, lamf
    double precision, parameter :: fac = 0.25d0 ! due to Gauss quadrature

    dxinv = 1.d0/dx

    Qflo = lo
    Qfhi(1:2) = hi(1:2)
    Qfhi(3) = hi(3)+1

    allocate(   Qf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3),QFVAR))
    allocate(dvdxf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3),2))
    allocate(dvdyf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3),2))
    allocate(  muf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3)))
    allocate(  xif(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3)))
    allocate( lamf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3)))
    allocate(Ddiaf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3),NSPEC))

    do n=1,QFVAR
       call cc2zface_3d(lo,hi,Q(:,:,:,n),qlo,qhi,Qf(:,:,:,n),Qflo,Qfhi)
    end do
    do n=1,2
       call cc2zface_3d(lo,hi,dveldx(:,:,:,n),qlo,qhi,dvdxf(:,:,:,n),Qflo,Qfhi)
    end do
    do n=1,2
       call cc2zface_3d(lo,hi,dveldy(:,:,:,n),qlo,qhi,dvdyf(:,:,:,n),Qflo,Qfhi)
    end do
    call cc2zface_3d(lo,hi, mu,qlo,qhi, muf,Qflo,Qfhi)
    call cc2zface_3d(lo,hi, xi,qlo,qhi, xif,Qflo,Qfhi)
    call cc2zface_3d(lo,hi,lam,qlo,qhi,lamf,Qflo,Qfhi)    
    do n=1,nspec
       call cc2zface_3d(lo,hi,Ddia(:,:,:,n),qlo,qhi,Ddiaf(:,:,:,n),Qflo,Qfhi)
    end do

    call comp_diff_flux_z(Qflo, Qfhi, f, flo, fhi, &
         Qf, muf, xif, lamf, Ddiaf, dvdxf, dvdyf, Qflo, Qfhi, &
         Q, qlo, qhi, dxinv, fac, domlo, domhi)

    deallocate(Qf, dvdxf, dvdyf, Ddiaf, muf, xif, lamf)

  end subroutine diff_z


  subroutine comp_diff_flux_z(lo, hi, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dveldx, dveldy, Qflo, Qfhi, &
       Q , Qlo, Qhi, dxinv, fac, domlo, domhi)

    use meth_params_module
    use derivative_stencil_module, only : FD4

    integer, intent(in) :: lo(3),hi(3),flo(3),fhi(3),Qflo(3),Qfhi(3),Qlo(3),Qhi(3),domlo(3),domhi(3)
    double precision, intent(in) :: dxinv(3), fac
    double precision, intent(inout)::   flx( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision, intent(in   )::    Qf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3),QFVAR)
    double precision, intent(in   )::    mu(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3))
    double precision, intent(in   )::    xi(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3))
    double precision, intent(in   )::   lam(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3))
    double precision, intent(in   )::  Ddia(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3),NSPEC)
    double precision, intent(in   )::dveldx(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3),2)
    double precision, intent(in   )::dveldy(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),Qflo(3):Qfhi(3),2)
    double precision, intent(in   )::     Q( Qlo(1): Qhi(1), Qlo(2): Qhi(2), Qlo(3): Qhi(3),QFVAR)

    integer :: i, j, k, n, UYN, QYN, QXN, QHN
    double precision :: tauzz, tauxz, tauyz
    double precision :: dudx, dudz, dvdy, dvdz, dwdx, dwdy, dwdz, divu
    double precision :: dTdz, dXdz, Vd
    double precision :: ek, rhovn
    double precision, dimension(lo(1):hi(1)) :: dlnpdz, Vc

    do k      =lo(3),hi(3)
       do j   =lo(2),hi(2)

          do i=lo(1),hi(1)
             ! viscous stress
             dudz = dxinv(3)*(FD4(-2)*Q(i,j,k-2,QU) + FD4(-1)*Q(i,j,k-1,QU) &
                  + FD4(0)*Q(i,j,k,QU) + FD4(1)*Q(i,j,k+1,QU))
             dvdz = dxinv(3)*(FD4(-2)*Q(i,j,k-2,QV) + FD4(-1)*Q(i,j,k-1,QV) &
                  + FD4(0)*Q(i,j,k,QV) + FD4(1)*Q(i,j,k+1,QV))
             dwdz = dxinv(3)*(FD4(-2)*Q(i,j,k-2,QW) + FD4(-1)*Q(i,j,k-1,QW) &
                  + FD4(0)*Q(i,j,k,QW) + FD4(1)*Q(i,j,k+1,QW))
             dudx = dxinv(1)*dveldx(i,j,k,1)
             dwdx = dxinv(1)*dveldx(i,j,k,2)
             dvdy = dxinv(2)*dveldy(i,j,k,1)
             dwdy = dxinv(2)*dveldy(i,j,k,2)
             divu = dudx + dvdy + dwdz
             tauxz = mu(i,j,k)*(dudz+dwdx)
             tauyz = mu(i,j,k)*(dvdz+dwdy)
             tauzz = mu(i,j,k)*(2.d0*dwdz-twoThirds*divu) + xi(i,j,k)*divu
             flx(i,j,k,UMX)   = flx(i,j,k,UMX)   - fac*tauxz
             flx(i,j,k,UMY)   = flx(i,j,k,UMY)   - fac*tauyz
             flx(i,j,k,UMZ)   = flx(i,j,k,UMZ)   - fac*tauzz
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac* &
                  (tauxz*Qf(i,j,k,QU) + tauyz*Qf(i,j,k,QV) + tauzz*Qf(i,j,k,QW))

             !thermal conduction
             dTdz = dxinv(3) * (FD4(-2)*Q(i,j,k-2,QTEMP) + FD4(-1)*Q(i,j,k-1,QTEMP) &
                  + FD4(0)*Q(i,j,k,QTEMP) + FD4(1)*Q(i,j,k+1,QTEMP))
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac*lam(i,j,k)*dTdz

             ! compute dpdz
             dlnpdz(i) = dxinv(3) * (FD4(-2)*Q(i,j,k-2,QPRES) + FD4(-1)*Q(i,j,k-1,QPRES) &
                  + FD4(0)*Q(i,j,k,QPRES) + FD4(1)*Q(i,j,k+1,QPRES)) / Qf(i,j,k,QPRES)
             Vc(i) = 0.d0
          end do

          do n=1,NSPEC
             UYN = UFS+n-1
             QYN = QFY+n-1
             QXN = QFX+n-1
             QHN = QFH+n-1
             do i = lo(1), hi(1)
                dXdz = dxinv(3) * (FD4(-2)*Q(i,j,k-2,QXN) + FD4(-1)*Q(i,j,k-1,QXN) &
                     + FD4(0)*Q(i,j,k,QXN) + FD4(1)*Q(i,j,k+1,QXN))
                Vd = -Ddia(i,j,k,n)*(dXdz + (Qf(i,j,k,QXN)-Qf(i,j,k,QYN))*dlnpdz(i))
                
                flx(i,j,k,UYN) = flx(i,j,k,UYN) + fac*Vd
                Vc(i) = Vc(i) + Vd
                flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) + fac*Vd*Qf(i,j,k,QHN)
             end do
          end do
    
          do n=1,NSPEC
             UYN = UFS+n-1
             QYN = QFY+n-1
             QHN = QFH+n-1
             do i = lo(1), hi(1)
                flx(i,j,k,UYN )  = flx(i,j,k,UYN  ) - (fac*Qf(i,j,k,QYN)*Vc(i))
                flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - (fac*Qf(i,j,k,QYN)*Vc(i))*Qf(i,j,k,QHN)
             end do
          end do
       end do
    end do

    if (.not. do_weno) then
       ! compute hyperbolic flux
       do k      =lo(3),hi(3)
          do j   =lo(2),hi(2)
             do i=lo(1),hi(1)
             rhovn = fac*Qf(i,j,k,QRHO)*Qf(i,j,k,QW)
             flx(i,j,k,URHO) = flx(i,j,k,URHO) + rhovn
             flx(i,j,k,UMX ) = flx(i,j,k,UMX ) + rhovn*Qf(i,j,k,QU)
             flx(i,j,k,UMY ) = flx(i,j,k,UMY ) + rhovn*Qf(i,j,k,QV)
             flx(i,j,k,UMZ ) = flx(i,j,k,UMZ ) + rhovn*Qf(i,j,k,QW) + fac*Qf(i,j,k,QPRES)

             ek = 0.5d0*(Qf(i,j,k,QU)**2+Qf(i,j,k,QV)**2+Qf(i,j,k,QW)**2)
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) + rhovn*ek
             do n=1,NSPEC
                flx(i,j,k,UEDEN)   = flx(i,j,k,UEDEN)   + (rhovn*Qf(i,j,k,QFY+n-1))*Qf(i,j,k,QFH+n-1)
                flx(i,j,k,UFS+n-1) = flx(i,j,k,UFS+n-1) + (rhovn*Qf(i,j,k,QFY+n-1))
             end do                
             end do
          end do
       end do
    end if
    
  end subroutine comp_diff_flux_z

end module difterm_module
