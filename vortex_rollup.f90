program vortex_rollup 
        implicit none
        integer, parameter :: sp =kind(1.0), dp=kind(1.0d0)
        integer, parameter :: wp = dp

        integer, parameter :: N=400
        real(kind=wp), parameter :: pi=3.14159265358979323_wp
        real(kind=wp), parameter :: dt=0.1_wp, del=0.5_wp
        real(kind=wp), dimension(N) :: x,y,xnew,ynew,xtemp,ytemp
        real(kind=wp), dimension(N) :: k1x,k2x,k3x,k4x, k1y,k2y,k3y,k4y
        integer :: numsteps=4/dt
        integer :: i

        do i=1,N
                x(i)=real(i-1,wp)/real(N,wp)  +0.01*sin(2*pi*real(i-1,wp)/real(n,wp))
                y(i)=-0.01*sin(2*pi*real(i-1,wp)/real(n,wp))
        end do
        do i=1,numsteps
                call rk4(x,y,N)
        end do
        do i=1,N
               print *, x(i), y(i)
        end do
contains 
        subroutine rhs(xnew,x,ynew,y,N)
                integer, intent(in) :: N
                real(kind=wp), intent(in), dimension(N) :: x,y
                real(kind=wp), intent(out),dimension(N) :: xnew,ynew
                integer :: j,k
                xnew=0_wp
                ynew=0_wp
                do j=1,N
                        do k=1,N
                             if(k /= j) then
                                  xnew(j)=xnew(j)+sinh(2*pi*(y(j)-y(k)))/(cosh(2*pi*(y(j)-y(k)))-cos(2*pi*(x(j)-x(k)))+del*del)
                                  ynew(j)=ynew(j)+sin(2*pi*(y(j)-y(k)))/(cosh(2*pi*(y(j)-y(k)))-cos(2*pi*(x(j)-x(k)))+del*del)
                              endif 
                        end do
                end do
                xnew=-xnew/real(2*N,wp)
                ynew=ynew/real(2*N,wp)
        end subroutine rhs
        subroutine rk4(x,y,N)
                integer, intent(in) :: N
                real(kind=wp),intent(inout), dimension(N):: x,y

                real(kind=wp), dimension(N) :: xnew,ynew,xtemp,ytemp
                real(kind=wp), dimension(N) :: k1x,k2x,k3x,k4x, k1y,k2y,k3y,k4y
                call rhs(k1x,x,k1y,y,N)
                k1x=k1x*dt
                k1y=k1y*dt
                xtemp=x+0.5_wp*k1x
                ytemp=y+0.5_wp*k1y
                call rhs(k2x,xtemp,k2y,ytemp,N)
                k2x=k2x*dt
                k2y=k2y*dt
                xtemp=x+0.5_wp*k2x
                ytemp=y+0.5_wp*k2y
                call rhs(k3x,xtemp,k3y,ytemp,N)
                k3x=k3x*dt
                k3y=k3y*dt
                xtemp=x+k3x
                ytemp=y+k3y
                call rhs(k4x,xtemp,k4y,ytemp,N) 

                xnew=x+(k1x+2.0_wp*k2x+2.0_wp*k3x+k4x)/real(6,wp)
                ynew=y+(k1y+2.0_wp*k2y+2.0_wp*k3y+k4y)/real(6,wp)
                x=xnew
                y=ynew
        end subroutine rk4
end program vortex_rollup

