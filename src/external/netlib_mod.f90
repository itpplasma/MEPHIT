module netlib_mod
  implicit none
  private
  public :: zeroin
contains

!  To get d1mach, mail netlib                                           
!       send d1mach from core                                           
      double precision function zeroin(ax,bx,f,tol) 
      double precision ax,bx,f,tol 
!                                                                       
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!                                                                       
!  input..                                                              
!                                                                       
!  ax     left endpoint of initial interval                             
!  bx     right endpoint of initial interval                            
!  f      function subprogram which evaluates f(x) for any x in         
!         the interval  ax,bx                                           
!  tol    desired length of the interval of uncertainty of the          
!         final result (.ge.0.)                                         
!                                                                       
!  output..                                                             
!                                                                       
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx    
!                                                                       
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not      
!  satisfied.   zeroin  returns a zero  x  in the given interval        
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.                                 
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).        
!                                                                       
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s 
      double precision  dabs, d1mach 
   10 eps = d1mach(4) 
      tol1 = eps+1.0d0 
!                                                                       
      a=ax 
      b=bx 
      fa=f(a) 
      fb=f(b) 
!     check that f(ax) and f(bx) have different signs                   
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20 
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20 
         write(6,2500) 
 2500    format(1x,'f(ax) and f(bx) do not have different signs,',      &
     &             ' zeroin is aborting')                               
         return 
   20 c=a 
      fc=fa 
      d=b-a 
      e=d 
   30 if (dabs(fc).ge.dabs(fb)) go to 40 
      a=b 
      b=c 
      c=a 
      fa=fb 
      fb=fc 
      fc=fa 
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol 
      xm = 0.5d0*(c-b) 
      if ((dabs(xm).le.tol1).or.(fb.eq.0.0d0)) go to 150 
!                                                                       
! see if a bisection is forced                                          
!                                                                       
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50 
      d=xm 
      e=d 
      go to 110 
   50 s=fb/fa 
      if (a.ne.c) go to 60 
!                                                                       
! linear interpolation                                                  
!                                                                       
      p=2.0d0*xm*s 
      q=1.0d0-s 
      go to 70 
!                                                                       
! inverse quadratic interpolation                                       
!                                                                       
   60 q=fa/fc 
      r=fb/fc 
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0)) 
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0) 
   70 if (p.le.0.0d0) go to 80 
      q=-q 
      go to 90 
   80 p=-p 
   90 s=e 
      e=d 
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge.            &
     &dabs(0.5d0*s*q))) go to 100                                       
      d=p/q 
      go to 110 
  100 d=xm 
      e=d 
  110 a=b 
      fa=fb 
      if (dabs(d).le.tol1) go to 120 
      b=b+d 
      go to 140 
  120 if (xm.le.0.0d0) go to 130 
      b=b+tol1 
      go to 140 
  130 b=b-tol1 
  140 fb=f(b) 
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20 
      go to 30 
  150 zeroin=b 
      return 
      END                                           

end module netlib_mod
