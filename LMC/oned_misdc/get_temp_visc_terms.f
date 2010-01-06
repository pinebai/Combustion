      subroutine get_temp_visc_terms(scal,beta,visc,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1)
      real*8 dx, time
      
      integer i
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 dxsqinv
      
c     Compute Div(lambda.Grad(T)) + rho.D.Grad(Hi).Grad(Yi)
C CEG:: the first thing rhoDgradHgradY() does is set grow cells
C      call set_bc_grow_s(scal,dx,time)
      call rhoDgradHgradY(scal,beta,visc,dx,time)

      open(UNIT=11, FILE='rdghgy.dat', STATUS = 'REPLACE')
      do i = 0,nx-1
         write(11,1006) i,1.d-1*visc(i)
      enddo
 1006 FORMAT((I5,1X),(E22.15,1X))    
      stop       



c     Add Div( lambda Grad(T) )
      dxsqinv = 1.d0/(dx*dx)
      do i = 0,nx-1
         if (coef_avg_harm.eq.1) then
            beta_lo = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i-1,Temp))
            beta_hi = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i+1,Temp))
         else
            beta_lo = 0.5*(beta(i,Temp) + beta(i-1,Temp))
            beta_hi = 0.5*(beta(i,Temp) + beta(i+1,Temp))
         endif
         
         flux_hi = beta_hi*(scal(i+1,Temp) - scal(i  ,Temp)) 
         flux_lo = beta_lo*(scal(i  ,Temp) - scal(i-1,Temp)) 
         visc(i) =  visc(i) + (flux_hi - flux_lo) * dxsqinv
      enddo
      end

      subroutine rhoDgradHgradY(scal,beta,visc,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1)
      real*8 dx,time
      
      integer i,n,is,IWRK
      real*8 beta_lo,beta_hi
      real*8 rdgydgh_lo,rdgydgh_hi
      real*8 dxsqinv,RWRK,rho,dv
      real*8 hi(maxspec,-1:nx)
      real*8 Y(maxspec,-1:nx)

C FIXME
      real*8 beta_edge(-1:nx  ,Nspec)      
      integer j

c     Compute rhoD Grad(Yi).Grad(hi) terms

      dxsqinv = 1.d0/(dx*dx)
c     Get Hi, Yi at cell centers
      call set_bc_grow_s(scal,dx,time)
      do i = -1,nx
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         call CKHMS(scal(i,Temp),IWRK,RWRK,hi(1,i))
         do n=1,Nspec
            Y(n,i) = scal(i,FirstSpec+n-1)/rho
C CEG:: changing tobe consistent with LMC
C            Y(n,i) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo
      enddo

C debugging FIXME
C$$$ 1007 FORMAT((I5,1X),6(E22.15,1X))      
C$$$         open(UNIT=11, FILE='h.dat', STATUS = 'REPLACE')
C$$$         write(11,*)'# 256 7' 
C$$$         do j=0,nx-1
C$$$            write(11,1007) j,(1.d-4*hi(n,j),n=1,Nspec)
C$$$         enddo
C$$$         close(11)
C$$$         stop
CCCCCCCCCCCCC


c     Compute differences
      do i = 0,nx-1
         dv = 0.d0
         do n=1,Nspec
            is = FirstSpec + n - 1
            if (coef_avg_harm.eq.1) then
               beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
               beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
            else
               beta_lo = 0.5*(beta(i,is) + beta(i-1,is))
               beta_hi = 0.5*(beta(i,is) + beta(i+1,is))
            endif
            
C CEG:: debugging remove me!!
            beta_edge(i,n) = beta_lo
            beta_edge(i+1,n) = beta_hi
CCCCC
            rdgydgh_lo = beta_lo*(Y(n,i)-Y(n,i-1))*(hi(n,i)-hi(n,i-1))
            rdgydgh_hi = beta_hi*(Y(n,i+1)-Y(n,i))*(hi(n,i+1)-hi(n,i))
            
            dv = dv + (rdgydgh_hi + rdgydgh_lo)*0.5d0
         enddo
         
         visc(i) = dv*dxsqinv
        end do

      open(UNIT=11, FILE='rdghgy.dat', STATUS = 'REPLACE')
      write(11,*)
      do i = 0,nx-1
         write(11,1006) i,1.d-1*visc(i)
      enddo
 1006 FORMAT((I5,1X),(E22.15,1X))    
      stop       

      end
