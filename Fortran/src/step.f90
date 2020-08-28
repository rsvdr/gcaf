!---------------------------------------------------------------------!
! GCAF= Guiding Center Analytic Field code by R. Saavedra             !
! Email of the author: saavestro@gmail.com                            !
! For more information please visit https://github.com/Saavestro/gcaf !
!---------------------------------------------------------------------!
!---------------------------------------------------------!
subroutine onestep(k)
!---------------------------------------------!
! Fourth order Runge-Kutta algorithm du/dt= e !
!---------------------------------------------!
!=========================================================!
	use shared,only: &
		pi,&    !pi number
		cs,&	!test type ions/elec
		pol,&   !poloidal flux
		tht,&   !theta
		zet,&   !zeta
		rho,&   !parallel gyroradius
		al,&    !alfa perturbation
		dadp,&  !al pol derivative
		dadt,&  !al tht derivative
		dadz,&  !al zet derivative
		dptdp,& !pt pol derivative
		dptdt,& !pt tht derivative
		dptdz,& !pt zet derivative
		dadtim,& !pt zet derivative
		dt,&    !time step
		qf,&    !q function
		bf,&    !b field
		dbdp,&  !b pol derivative
		dbdt,&  !b tht derivative
		dbdz,&  !b zet derivative
		cg,&    !poloidal current
		cgp,&   !cg pol derivative
		ci,&    !toroidal current
		cip,&   !ci pol derivative
		mu,&    !b moment
		polw,&  !on wall poloidal flux
		dto     !initial time step
	implicit none
	integer(8):: k
	integer(8):: i
	integer(8):: lcl
	integer(8):: ndum
	real(8):: e(4)
	real(8):: u0(4)
	real(8):: u1(4)
	real(8):: f0(4)
	real(8):: f1(4)
	real(8):: f2(4)
	real(8):: f3(4)
	real(8):: xdum
	real(8):: ydum
	real(8):: xdot
	real(8):: ydot
	real(8):: rsq
	real(8):: fac1
	real(8):: fac2
	real(8):: fac3
	real(8):: det
	real(8):: rbb
	real(8):: dedb
!=========================================================!
	lcl= .6d0*(1.d0 +sign(1.d0,pol(k) -5.d-2*polw))
	!lcl= 1
	dt(k)= lcl*dt(k) +(1 -lcl)*dto
	xdum= sqrt(pol(k))*cos(tht(k))
	ydum= sqrt(pol(k))*sin(tht(k))
	u1(1)= lcl*pol(k)+ (1- lcl)*xdum
	u1(2)= lcl*tht(k)+ (1- lcl)*ydum
	u1(3)= zet(k)
	u1(4)= rho(k)
	u0(1)= u1(1)
	u0(2)= u1(2)
	u0(3)= u1(3)
	u0(4)= u1(4)

	do i=1,4
		call field(k)
		rbb= u1(4)*bf(k)*bf(k)
		dedb= mu(k) +u1(4)*u1(4)*bf(k)
		det= cg(k)*qf(k) +ci(k) &
			+(cs*u1(4) +al(k))*(cg(k)*cip(k) -ci(k)*cgp(k))
		fac1= 1.d0 -cgp(k)*(cs*u1(4) +al(k)) -cg(k)*dadp(k)
		fac2= qf(k) +cip(k)*(cs*u1(4) +al(k)) +ci(k)*dadp(k)
		fac3= cs*ci(k)*dadz(k)- cs*cg(k)*dadt(k)
		e(1)= -cs*cg(k)*dedb*dbdt(k)/det &
			+cs*ci(k)*dedb*dbdz(k)/det &
			+cg(k)*rbb*dadt(k)/det &
			-ci(k)*rbb*dadz(k)/det &
			-cs*cg(k)*dptdt(k)/det &
			+cs*ci(k)*dptdz(k)/det
		e(2)= cs*cg(k)*dedb*dbdp(k)/det &
			+rbb*fac1/det &
			+cs*cg(k)*dptdp(k)/det
		e(3)= -cs*ci(k)*dedb*dbdp(k)/det &
			+rbb*fac2/det &
			-cs*ci(k)*dptdp(k)/det
		e(4)= -fac1*dedb*dbdt(k)/det &
			-fac2*dedb*dbdz(k)/det &
			+fac3*dedb*dbdp(k)/det &
			-fac1*dptdt(k)/det &
			-fac2*dptdz(k)/det &
			+fac3*dptdp(k)/det &
			-dadtim(k)
		rsq= u1(1)*u1(1)+ u1(2)*u1(2)
		xdot= .5d0*u1(1)*e(1)/rsq -u1(2)*e(2)
		ydot= .5d0*u1(2)*e(1)/rsq +u1(1)*e(2)
		e(1)= lcl*e(1) +(1 -lcl)*xdot
		e(2)= lcl*e(2) +(1 -lcl)*ydot

		go to(11,12,13,14),i
		11 continue
		f0(1:4)= e(1:4)
		u1(1:4)= u0(1:4)+ .5d0*dt(k)*f0(1:4)
		go to 40
		12 continue
		f1(1:4)= e(1:4)
		u1(1:4)= u0(1:4)+ .5d0*dt(k)*f1(1:4)
		go to 40
		13 continue
		f2(1:4)= e(1:4)
		u1(1:4)= u0(1:4)+ dt(k)*e(1:4)
		go to 40
		14 continue
		f3(1:4)= e(1:4)
		u1(1:4)= u0(1:4) +(dt(k)/6.d0)*(f0(1:4) +2.d0*f1(1:4) &
			+2.d0*f2(1:4) +f3(1:4))
		go to 40
		40 continue
	enddo

	ndum= .6d0*(1.d0 -sign(1.d0,u1(1)))
	pol(k)= lcl*u1(1) +(1- lcl)*(u1(1)*u1(1) +u1(2)*u1(2))
	tht(k)= lcl*u1(2) +(1- lcl)*(atan(u1(2)/u1(1)) +ndum*pi)
	zet(k)= u1(3)
	rho(k)= u1(4)
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine stepfix(k)
!-----------------------------------------!
! Tune the time step along the simulation !
!-----------------------------------------!
!=========================================================!
	use shared,only: &
		dt,&     !time step
		deler,&  !input step fix en delta
		totenn,& !total en normalized
		dto      !initial time step
	implicit none
	integer(8):: k
	integer(8):: ndum
	real(8):: deler2
	real(8):: dum
!=========================================================!
	deler2= 1.d-1*deler
	ndum= .6d0*(1 +sign(1.d0,totenn(k)/deler -1.d0)) !1 if totenn > dele
	dt(k)= dt(k)*(1 -ndum) +.8d0*dt(k)*ndum
	ndum= .6d0*(1 -sign(1.d0,totenn(k)/deler2 -1.d0)) !1 if totenn < dele2
	dt(k)= dt(k)*(1 -ndum) +1.2d0*dt(k)*ndum
	dum= 2.5d-1*dto
	dt(k)= min(dum,dt(k))
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine update(k,nskip,mskip,nend)
!-----------------------------!
! Update the integration step !
!-----------------------------!
!=========================================================!
	use shared,only: &
		keverg,&  !keV to erg
		ndisp,&   !input number of disp points
		ndist,&   !input number of dist points
		rmin,&    !minor radius
		en,&	  !energy
		ptch,&    !pitch
		pol,&     !poloidal flux
		tht,&     !theta
		zet,&     !zeta
		rho,&     !parallel gyroradius
		toten,&   !total en
		totenn,&  !total en normalized
		time,&    !time
		dt,&      !time step
		dto,&     !initial time step
		pz,&      !zeta conjugate momenta
		cg,&      !poloidal current
		al,&      !alpha perturbation
		bf,&      !b field
		pt,&      !electrostatic potential
		mu,&      !b moment
		trun,&    !total simulation time
		enon,&    !normalized initial en
		enaxo,&   !on axis en
		polw,&    !on wall pol
		leno,&    !characteristic length
		omego,&   !characteristic time freq
		dft,&     !time for disp
		rprv,&    !previous rad
		flux,&    !particle flux
		rone,&    !rad one part
		rsqone,&  !rad sq one part
		rcuone,&  !rad cu one part
		rfoone,&  !rad fo one part
		enone,&   !en one part
		ensqone,& !en sq one part
		encuone,& !en cu one part
		enfoone,& !en fo one part
		pzone,&   !pz one part
		pzsqone,& !pz sq one part
		pzcuone,& !pz cu one part
		pzfoone,& !pz fo one part
		distt,&   !time dist
		distr,&   !rad dist
		diste,&   !en dist
		distp     !ptch dist
	implicit none
	integer(8):: k
	integer(8):: nskip
	integer(8):: mskip
	integer(8):: nend
	real(8):: dtdum    !dt dummy
	real(8):: rdum     !r dummy
	real(8):: raddum   !r/rmin dummy
	real(8):: flxpln   !flux plane
	real(8):: prvdum   !rprv/rmin dummy
!=========================================================!
	call field(k)
	toten(k)= .5d0*bf(k)*bf(k)*rho(k)*rho(k) +mu(k)*bf(k) +pt(k)
	totenn(k)= abs((toten(k) -en(k))/enon)
	en(k)= toten(k) -pt(k)
	ptch(k)= rho(k)*bf(k)/sqrt(2*en(k))
	time(k)= time(k) +dt(k)
	pz(k)= cg(k)*(rho(k) +al(k)) -pol(k)

	flxpln= 0.52d0
	raddum= leno*sqrt(2*pol(k))/rmin
	prvdum= rprv(k)/rmin
	if(raddum.gt.flxpln.and.prvdum.lt.flxpln)then
		flux= flux +1
	endif
	if(raddum.lt.flxpln.and.prvdum.gt.flxpln)then
		flux= flux -1
	endif
	rprv(k)= leno*sqrt(2*pol(k))

	if(time(k).gt.nskip*trun/ndisp) then
		rdum= leno*sqrt(2*pol(k))
		dft(nskip)= time(k)
		rone(nskip)= rdum
		rsqone(nskip)= rdum*rdum
		rcuone(nskip)= rdum*rdum*rdum
		rfoone(nskip)= rdum*rdum*rdum*rdum
		enone(nskip)= en(k)
		ensqone(nskip)= en(k)*en(k)
		encuone(nskip)= en(k)*en(k)*en(k)
		enfoone(nskip)= en(k)*en(k)*en(k)*en(k)
		pzone(nskip)= pz(k)
		pzsqone(nskip)= pz(k)*pz(k)
		pzcuone(nskip)= pz(k)*pz(k)*pz(k)
		pzfoone(nskip)= pz(k)*pz(k)*pz(k)*pz(k)
		nskip= nskip +1
	endif

	if(time(k).gt.mskip*trun/ndist) then
		distt(mskip)= time(k)/omego
		distr(k,mskip)= leno*sqrt(2*pol(k))
		diste(k,mskip)= enaxo*en(k)/keverg
		distp(k,mskip)= ptch(k)
		mskip= mskip +1
	endif

	if(time(k).gt.trun) then
		nend= 1
		return
	endif
	dtdum= 1.d-6*dto
	if(dt(k).lt.dtdum) then
		write(0,*) 'small step dt'
		nend= 1
	endif
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine update_ensemble
!---------------------------!
! Update ensemble variables !
!---------------------------!
!=========================================================!
	use shared,only: &
		rav,&     !rad av
		rsqav,&   !rad sq av
		rcuav,&   !rad cu av
		rfoav,&   !rad fo av
		enav,&    !en av
		ensqav,&  !en sq av
		encuav,&  !en cu av
		enfoav,&  !en fo av
		pzav,&    !pz av
		pzsqav,&  !pz sq av
		pzcuav,&  !pz cu av
		pzfoav,&  !pz fo av
		rone,&    !rad one part
		rsqone,&  !rad sq one part
		rcuone,&  !rad cu one part
		rfoone,&  !rad fo one part
		enone,&   !en one part
		ensqone,& !en sq one part
		encuone,& !en cu one part
		enfoone,& !en fo one part
		pzone,&   !pz one part
		pzsqone,& !pz sq one part
		pzcuone,& !pz cu one part
		pzfoone,& !pz fo one part
		ndisp	  !input number of disp points
	implicit none
!=========================================================!
	rav(1:ndisp)= rav(1:ndisp)+ rone(1:ndisp)
	rsqav(1:ndisp)= rsqav(1:ndisp)+ rsqone(1:ndisp)
	rcuav(1:ndisp)= rcuav(1:ndisp)+ rcuone(1:ndisp)
	rfoav(1:ndisp)= rfoav(1:ndisp)+ rfoone(1:ndisp)
	enav(1:ndisp)= enav(1:ndisp)+ enone(1:ndisp)
	ensqav(1:ndisp)= ensqav(1:ndisp)+ ensqone(1:ndisp)
	encuav(1:ndisp)= encuav(1:ndisp)+ encuone(1:ndisp)
	enfoav(1:ndisp)= enfoav(1:ndisp)+ enfoone(1:ndisp)
	pzav(1:ndisp)= pzav(1:ndisp)+ pzone(1:ndisp)
	pzsqav(1:ndisp)= pzsqav(1:ndisp)+ pzsqone(1:ndisp)
	pzcuav(1:ndisp)= pzcuav(1:ndisp)+ pzcuone(1:ndisp)
	pzfoav(1:ndisp)= pzfoav(1:ndisp)+ pzfoone(1:ndisp)
	return
end subroutine
!---------------------------------------------------------!
