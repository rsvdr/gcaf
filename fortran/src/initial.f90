!---------------------------------------------------------------------!
! GCAF= Guiding Center Analytic Field code by R. Saavedra             !
! For more information please visit https://git.rsaavedra.xyz/?p=gcaf.git !
!---------------------------------------------------------------------!
!---------------------------------------------------------!
subroutine inivar
!----------------------!
! Initialize variables !
!----------------------!
!=========================================================!
	use shared,only: &
		pi,&       !pi number
		keverg,&   !keV to erg
		btkg,&     !10 kG
		rq,&       !electron charge
		rc,&       !speed of light
		mp,&       !proton mass
		cs,&       !test type ions/elec
		mt,&       !test mass in protons mass
		zt,&       !test charge in protons charge
		ntype,&    !simulation type traj/poin/ense
		ntran,&    !number of toroidal transits
		nsdel,&    !input step delta
		npart,&    !number of particles
		nengo,&    !initial en in keV
		nrado,&    !initial radius from b axis
		npto,&     !electrostatic potential in keV
		ptch,&     !pitch angle
		leno,&     !characteristic length
		omego,&    !characteristic time frequency
		dto,&      !initial time step
		pto,&      !pt magnitude
		zplo,&     !poin zeta plane
		bm,&       !on axis b field magnitude
		rmaj,&     !major radius
		rmin,&     !minor radius
		enon,&     !normalized initial en
		enaxo,&    !on axis particle en
		ptchgrad,& !pitch in grades
		thtgrad,&  !tht in grades
		zetgrad,&  !zet in grades
		polw,&     !on wall poloidal flux
		trun,&     !total simulation time
		ptcho,&    !initial pitch angle in cos
		eno,&      !initial en in keV
		polo,&     !initial poloidal flux
		theto,&    !initial theta
		zeto,&     !initial zeta
		nfrq,&     !input perturbation frequency
		prtfrq,&   !perturbation frequency
		flux       !particle flux
	implicit none
	integer(8):: mdum !dummy zet
	real(8):: rdum    !dummy rad variable
	real(8):: tran 	  !number of transits
	real(8):: go      !poloidal current
	real(8):: bt      !toroidal field
!=========================================================!

	if(ntype.eq.1) npart= 1
	leno= 1.d0*rmaj 
	polw= .5d0*rmin*rmin/leno/leno
	theto= thtgrad*pi/180
	zeto= zetgrad*pi/360
	rdum= 1.d0*nrado
	if(nrado.eq.0) rdum= rmin/2
	polo= .5d0*rdum*rdum/leno/leno
	ptcho= cos(ptchgrad*pi/180)

	zplo= 0*pi/8 +pi/1000
	mdum= zplo/2/pi
	zplo= zplo -2*pi*mdum
	if(zplo.lt.0.d0) zplo= zplo +2*pi

	if(cs.eq.1) omego= rq*btkg*zt*bm/mt/mp/rc
	if(cs.eq.-1) omego= 1836.1526d0*rq*btkg*zt*bm/mt/mp/rc
	if(cs.eq.1) enaxo= mp*mt*leno*leno*omego*omego
	if(cs.eq.-1) enaxo= mp*mt*leno*leno*omego*omego/1836.1526d0
	prtfrq= nfrq/omego
	pto= keverg*npto/enaxo

	bt= bm
	eno= 1.d0*nengo
	enon= keverg*eno/enaxo
	go= rmaj*bt/leno/bm
	tran= 2*pi*go/bt/sqrt(2*enon)
	trun= ntran*tran
	dto= tran/nsdel
	flux= 0
	call rcrd_log(36)
	call info
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine set(k)
!----------------------!
! Set vector variables !
!----------------------!
!=========================================================!
	use shared,only: &
		pi,&     !pi number
		keverg,& !keV to erg
		mp,&     !proton mass
		mt,&     !test mass in protons mass
		ntype,&  !simulation type traj/poin/ense
		npart,&  !number of particles
		nends,&  !en distribution mon/max
		nrado,&  !initial radius in cm
		ptch,&   !pitch angle
		en,&     !energy
		pol,&    !poloidal flux
		tht,&    !theta
		zet,&    !zeta
		mu,&     !b moment
		pz,&     !zeta conjugate momenta
		rho,&    !parallel gyroradius
		r1,&     !poin lower bound
		r2,&     !poin upper bound
		leno,&   !characteristic length
		dt,&     !time step
		dto,&    !initial time step
		time,&   !time variable
		al,&     !alfa perturbation
		bf,&     !b field magnitude
		pt,&     !electrostatic potential
		bm,&     !on axis b
		rmaj,&   !major radius
		rmin,&   !minor radius
		enon,&   !normalized initial en
		enaxo,&  !on axis particle en
		zplo,&   !poin zeta plane
		ptcho,&  !initial pitch angle in cos
		eno,&    !initial en in keV
		polo,&   !initial poloidal flux
		theto,&  !initial theta
		zeto,&   !initial zeta
		cg       !poloidal current
	implicit none
	integer(8):: k
	real(8):: ranx    !random function
	real(8):: norm_ab !gaussian function
	real(8):: vx      !vel x component
	real(8):: vy      !vel y component
	real(8):: vz      !vel z component
	real(8):: endum   !en dummy
	real(8):: rdum    !r dummy
	real(8):: vsq     !vel sq dummy
	real(8):: p1      !poin pol lower bound
	real(8):: p2      !poin pol upper bound
	real(8):: r1dum   !poin rad lower bound dummy
	real(8):: r2dum   !poin rad upper bound dummy
!=========================================================!
	!set pol,thet,zet
	if(ntype.eq.1) then
		p1= polo
		p2= polo
	endif
	if(ntype.eq.2) then
		p1= polo
		p2= 1.01d0*polo
	endif
	if(ntype.eq.3) then
		r1dum= r1*rmin
		r2dum= r2*rmin
		p1= .5d0*r1dum*r1dum/leno/leno
		p2= .5d0*r2dum*r2dum/leno/leno
		rdum= r1dum+ k*(r2dum- r1dum)/npart
		pol(k)= polo
		if(nrado.eq.-1) pol(k)= .5d0*rdum*rdum/leno/leno
	endif
	if(ntype.ne.3) pol(k)= polo
	tht(k)= ranx()*2*pi
	tht(k)= theto
	zet(k)= ranx()*2*pi
	zet(k)= zeto
	if(ntype.eq.1) then
		pol(k)= polo
		tht(k)= theto
		zet(k)= zeto
	endif

	endum= eno
	vx= sqrt(keverg*endum/mp/mt)*norm_ab(0.d0,1.d0)
	vy= sqrt(keverg*endum/mp/mt)*norm_ab(0.d0,1.d0)
	vz= sqrt(keverg*endum/mp/mt)*norm_ab(0.d0,1.d0)
	vsq= vx*vx +vy*vy +vz*vz
	if(nends.eq.1)then
		endum= mp*mt*vsq/2/keverg
	endif
	enon= keverg*endum/enaxo

	!set t,dt,en,rho,mu
	call field(k) !needs pol,tht,zet
	time(k)= 0.d0
	dt(k)= dto
	if(ntype.eq.1) ptch(k)= ptcho
	if(ntype.eq.2) ptch(k)= sign(1.d0,0.5d0 -ranx())*ranx()
	if(ntype.eq.3) ptch(k)= ptcho
	en(k)= enon
	rho(k)= ptch(k)*sqrt(2.d0*en(k))/bf(k)
	mu(k)= en(k)/bf(k)- .5d0*rho(k)*rho(k)*bf(k)
	pz(k)= cg(k)*(rho(k) +al(k)) -pol(k)
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine inivec
!-----------------!
! Allocate memory !
!-----------------!
	use shared
	implicit none
	integer(8):: ierr
	allocate(dt(idm),stat=ierr)
	allocate(time(idm),stat=ierr)
	allocate(en(idm),stat=ierr)
	allocate(ptch(idm),stat=ierr)
	allocate(dft(idm),stat=ierr)
	allocate(pol(idm),stat=ierr)
	allocate(tht(idm),stat=ierr)
	allocate(zet(idm),stat=ierr)
	allocate(rho(idm),stat=ierr)
	allocate(toten(idm),stat=ierr)
	allocate(totenn(idm),stat=ierr)
	allocate(pold(idm),stat=ierr)
	allocate(told(idm),stat=ierr)
	allocate(zold(idm),stat=ierr)
	allocate(distr(idm,idm),stat=ierr)
	allocate(diste(idm,idm),stat=ierr)
	allocate(distp(idm,idm),stat=ierr)
	allocate(rcrdtraj(idm,idm),stat=ierr)
	allocate(rcrdpoin(idm,idm),stat=ierr)
	allocate(distt(ndist),stat=ierr)
	allocate(bf(idm),stat=ierr)
	allocate(cg(idm),stat=ierr)
	allocate(cgp(idm),stat=ierr)
	allocate(ci(idm),stat=ierr)
	allocate(cip(idm),stat=ierr)
	allocate(qf(idm),stat=ierr)
	allocate(dbdp(idm),stat=ierr)
	allocate(dbdt(idm),stat=ierr)
	allocate(dbdz(idm),stat=ierr)
	allocate(pt(idm),stat=ierr)
	allocate(dptdp(idm),stat=ierr)
	allocate(dptdt(idm),stat=ierr)
	allocate(dptdz(idm),stat=ierr)
	allocate(al(idm),stat=ierr)
	allocate(dadp(idm),stat=ierr)
	allocate(dadt(idm),stat=ierr)
	allocate(dadz(idm),stat=ierr)
	allocate(dadtim(idm),stat=ierr)
	allocate(mu(idm),stat=ierr)
	allocate(pz(idm),stat=ierr)
	allocate(rprv(idm),stat=ierr)
	allocate(rav(idm),stat=ierr)
	allocate(rsqav(idm),stat=ierr)
	allocate(rcuav(idm),stat=ierr)
	allocate(rfoav(idm),stat=ierr)
	allocate(enav(idm),stat=ierr)
	allocate(ensqav(idm),stat=ierr)
	allocate(encuav(idm),stat=ierr)
	allocate(enfoav(idm),stat=ierr)
	allocate(pzav(idm),stat=ierr)
	allocate(pzsqav(idm),stat=ierr)
	allocate(pzcuav(idm),stat=ierr)
	allocate(pzfoav(idm),stat=ierr)
	allocate(rone(idm),stat=ierr)
	allocate(rsqone(idm),stat=ierr)
	allocate(rcuone(idm),stat=ierr)
	allocate(rfoone(idm),stat=ierr)
	allocate(enone(idm),stat=ierr)
	allocate(ensqone(idm),stat=ierr)
	allocate(encuone(idm),stat=ierr)
	allocate(enfoone(idm),stat=ierr)
	allocate(pzone(idm),stat=ierr)
	allocate(pzsqone(idm),stat=ierr)
	allocate(pzcuone(idm),stat=ierr)
	allocate(pzfoone(idm),stat=ierr)
	if(ierr.gt.0) then
		print *, 'INSUFFICIENT MEMORY'
		stop
	endif
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine finvec
!-------------------!
! Deallocate memory !
!-------------------!
	use shared
	implicit none
	deallocate(dt)
	deallocate(time)
	deallocate(en)
	deallocate(ptch)
	deallocate(dft)
	deallocate(pol)
	deallocate(tht)
	deallocate(zet)
	deallocate(rho)
	deallocate(toten)
	deallocate(totenn)
	deallocate(pold)
	deallocate(told)
	deallocate(zold)
	deallocate(distr)
	deallocate(diste)
	deallocate(distp)
	deallocate(trgttraj)
	deallocate(trgtpoin)
	deallocate(rcrdtraj)
	deallocate(rcrdpoin)
	deallocate(trgtcoef)
	deallocate(bf)
	deallocate(cg)
	deallocate(cgp)
	deallocate(ci)
	deallocate(cip)
	deallocate(qf)
	deallocate(dbdp)
	deallocate(dbdt)
	deallocate(dbdz)
	deallocate(pt)
	deallocate(dptdp)
	deallocate(dptdt)
	deallocate(dptdz)
	deallocate(al)
	deallocate(dadp)
	deallocate(dadt)
	deallocate(dadz)
	deallocate(dadtim)
	deallocate(mu)
	deallocate(pz)
	deallocate(rprv)
	deallocate(rav)
	deallocate(rsqav)
	deallocate(rcuav)
	deallocate(rfoav)
	deallocate(enav)
	deallocate(ensqav)
	deallocate(encuav)
	deallocate(enfoav)
	deallocate(pzav)
	deallocate(pzsqav)
	deallocate(pzcuav)
	deallocate(pzfoav)
	deallocate(rone)
	deallocate(rsqone)
	deallocate(rcuone)
	deallocate(rfoone)
	deallocate(enone)
	deallocate(ensqone)
	deallocate(encuone)
	deallocate(enfoone)
	deallocate(pzone)
	deallocate(pzsqone)
	deallocate(pzcuone)
	deallocate(pzfoone)
	return
end subroutine
!---------------------------------------------------------!
