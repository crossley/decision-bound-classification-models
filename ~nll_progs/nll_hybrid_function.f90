program nll_hybrid_function

use anorin_int

implicit none

real :: a1,a2,b,xc,noise1,noise2
real :: a1New,a2New,bNew,xcNew

integer :: i,j,l,grid_width ! The maximum width of the z-grid is 1000
real :: p, stepsize
real, dimension(1000) :: z

integer, dimension(1000) :: Aindices, Bindices

integer :: numpoints, icount, numA, numB ! Note that 10000 is the max number of data points

double precision, dimension(10000) :: likeA, likeB, negloglikeA, negloglikeB
double precision :: negloglike, gridsize

real :: datax, datay

! We start by reading in the values for the parameters that'll we'll use later

open(10, file = 'params.dat')

! X = (m,b,xc,noise1,noise2)

read(10,*) a1, b, xc, noise1, noise2

close(10)

a2 = sqrt(1-a1*a1)

! Now we compute a 1-dimensional z-grid.  For this we read in the size of the grid.

open(20, file = 'function_info.dat')

read(20,*) grid_width

stepsize = 1.0/(grid_width)
p = (0.5)*stepsize
z = 0

do i = 1,grid_width

	z(i) = anorin(p)
	p = p + stepsize

end do

gridsize = grid_width*grid_width

! Next we read in the number of points, as well as the number of A and B responses.

read(20,*) numpoints, numA, numB

! Next we read in the indices of the A and B reponses

read(20,*) Aindices
read(20,*) Bindices

! Now we compute the negative log likelihood for the pattern of responses

do l = 1,numpoints

	read(20,*) datax, datay
	
	! Next we transform the bounds
	
	a1New = (a1*noise1)
	a2New = (a2*noise2)
	bNew = a1*datax+a2*datay+b
	xcNew = (xc-datax)/noise1
	
	! Now we are ready to compute the integral.
	
	icount = 0 
		
	! We step through the z-grid and make our calculation.
	
	do i = 1, grid_width
	
		do j = 1, grid_width
			
			if ((a1New*z(i)+a2New*z(j)+bNew > 0).and.(z(i)<xcNew)) icount = icount + 2
			if ((a1New*z(i)+a2New*z(j)+bNew == 0).and.(z(i)<xcNew)) icount = icount + 1
			if ((a1New*z(i)+a2New*z(j)+bNew > 0).and.(z(i)==xcNew)) icount = icount + 1
			if ((a1New*z(i)+a2New*z(j)+bNew == 0).and.(z(i)==xcNew)) icount = icount + 1		
	
		end do
		
	end do

	likeB(l) = cap((1.0*icount)/(2*gridsize))
	likeA(l) = cap(1-likeB(l))
	
end do

! We compute the negative log likelihood for the data set

negloglikeB = -log(likeB)
negloglikeA = -log(likeA)

negloglike = 0

do i = 1, numB
	
	negloglike = negloglike + negloglikeB(Bindices(i))
	
end do

do i = 1, numA

	negloglike = negloglike + negloglikeA(Aindices(i))
	
end do
 
open(30, file = 'negloglike.dat')

write(30,110) negloglike

 close(20)
 close(30)

100 format(2(i3,x),f6.4)
110 format(f32.16)

contains

function cap(input)

double precision :: cap

double precision :: input

if (input < 0.0013) input = 0.0013

if (input > 0.9987) input = 0.9987

 cap = input

end function cap

end program nll_hybrid_function

