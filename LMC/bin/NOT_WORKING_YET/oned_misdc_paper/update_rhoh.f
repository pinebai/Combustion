      subroutine update_rhoh(scal_old,scal_new,aofs,alpha,
     &                       beta,dRhs,Rhs,dx,dt,be_cn_theta,time)
      implicit none
      include 'spec.h'
      real*8 scal_old(-1:nx  ,nscal)
      real*8 scal_new(-1:nx  ,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8    alpha(0 :nx-1)
      real*8     beta(-1:nx  ,nscal)
      real*8     dRhs(0 :nx-1)
      real*8      Rhs(0 :nx-1)
      real*8 dx,dt,be_cn_theta,time

      real*8  visc(0:nx-1)
      real*8  h_hi,h_lo,h_mid
      real*8  flux_lo,flux_hi
      real*8  dxsqinv
      real*8  beta_lo,beta_hi
      real*8  visc_term, RWRK
      integer i,n,is, IWRK

      real*8 cp,T,H,rho
      real*8 Y(Nspec)

      call get_rhoh_visc_terms(scal_old,beta,visc,dx,time)
      
      do i = 0,nx-1
         visc_term = dt*(1.d0 - be_cn_theta)*visc(i)
         
         scal_new(i,RhoH) = scal_old(i,RhoH) + dt*aofs(i,RhoH)
         Rhs(i) = dRhs(i) + scal_new(i,RhoH) + visc_term
         alpha(i) = scal_new(i,Density)
      enddo
      
      end
