!---------------------------------------------------------------------!
! GCAF= Guiding Center Analytic Field code by R. Saavedra             !
! Email of the author: saavestro@gmail.com                            !
! For more information please visit https://github.com/Saavestro/gcaf !
!---------------------------------------------------------------------!
!---------------------------------------------------------!
real(8) function xxproj(px,tx,zx)
	use shared,only:leno,rmaj
	implicit none
	real(8):: px,tx,zx,pdum,rdum
	pdum= 2.d0*leno*leno*px
	rdum= sqrt(pdum)
	xxproj= (1.d0*rmaj +rdum*cos(tx))*cos(zx)
	return
end function
!---------------------------------------------------------!


!---------------------------------------------------------!
real(8) function yyproj(px,tx,zx)
	use shared,only:leno,rmaj
	implicit none
	real(8):: px,tx,zx,pdum,rdum
	pdum= 2.d0*leno*leno*px
	rdum= sqrt(pdum)
	yyproj= (1.d0*rmaj +rdum*cos(tx))*sin(zx)
	return
end function
!---------------------------------------------------------!


!---------------------------------------------------------!
real(8) function zzproj(px,tx,zx)
	use shared,only:leno
	implicit none
	real(8):: px,tx,zx,pdum,rdum
	pdum= 2.d0*leno*leno*px
	rdum= sqrt(pdum)
	zzproj= rdum*sin(tx)
	return
end function
!---------------------------------------------------------!


!---------------------------------------------------------!
real(8) function bfield(pdum,tdum,zdum)
	use shared,only:rmaj,rmin,polw
	implicit none
	real(8):: ewdum,pdum,tdum,zdum
	ewdum= 1.d0*rmin/rmaj
	bfield= 1.d0 -ewdum*sqrt(pdum/polw)*cos(tdum)
	return
end function
!---------------------------------------------------------!


!---------------------------------------------------------!
real(8) function ranx()
	implicit none
	real(8):: r
	call random_number(r)
	ranx= r
	return
end function
!---------------------------------------------------------!


!---------------------------------------------------------!
function norm_ab(a,b)
!----------------------------------------------!
! Gaussian function with mean 'a' and standard !
! deviation 'b'                                !
!----------------------------------------------!
	implicit none
	real(8) :: a,b,r1,r2,x,norm_ab,ranx
	real(8), parameter :: pi = 3.14159265d0
	r1= ranx()
	r2= ranx()
	x= sqrt(-2.0d0*log(r1))*cos(2.0d0*pi*r2)
	norm_ab= a +b*x
end function
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine least_square(n,x,y,b,d)
!-------------------------------------------------!
! Linear least squares subroutine                 !
! The input data set is X(m), Y(m). The number of !
! data points is n (n must be > 2). The returned  !
! parameters are: a,b, coefficients of equation   !
! Y = a + b X, and d, standard deviation of fit.  !
!-------------------------------------------------!
	implicit none
	integer(4):: n
	real(8)::  a,b,d,x(0:n),y(0:n)  
	real(8)::  a1,a2,b0,b1,d1
	integer*4 m
	a1= 0
	a2= 0
	b0= 0
	b1= 0
	do m= 0,n-1
		a1= a1+ x(m)
		a2= a2+ x(m)*x(m)
		b0= b0+ y(m)
		b1= b1+ y(m)*x(m)
	end do
	a1= a1/n
	a2= a2/n
	b0= b0/n
	b1= b1/n
	d= a1*a1 -a2
	a= a1*b1 -a2*b0
	a= a/d
	b= a1*b0 -b1
	b= b/d
	!evaluation of standard deviation d (unbiased estimate) 
	d= 0
	do m= 0, n-1
		d1= y(m) -a -b*x(m)
		d= d +d1*d1
	end do
	d= sqrt(d/(n -2))
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine strip(long,string,n)
!---------------------------------------------------!
! strip long into string where there are n elements !
!---------------------------------------------------!
	use shared
	implicit none
	character(64):: string
	character(64):: long,str_old,text
	integer(4):: stringLen 
	integer(4):: last, actual
	integer(8):: n,i

	string= long(9:) !cut user input string

	!Strip string 
	stringLen = len (string)
	last = 1
	actual = 1

	do while (actual < stringLen)
		if (string(last:last) == ' ') then
		    actual = actual + 1
		    string(last:last) = string(actual:actual)
		    string(actual:actual) = ' '
		else
		    last = last + 1
		    if (actual < last) &
		        actual = last
		endif
	enddo

	n= len_trim(string)/3 !n= number of elements in strp
	allocate(prnt(n))
	read(long,*) text,prnt(1:n) !prnt= character vector
	str_old= prnt(1) 
	do i=2,size(prnt) !str= prnt concatenation
		string= trim(str_old)//trim(prnt(i))
		str_old= trim(string)
	enddo
	deallocate(prnt)
	return
end subroutine
!---------------------------------------------------------!
