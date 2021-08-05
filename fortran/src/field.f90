!---------------------------------------------------------------------!
! GCAF= Guiding Center Analytic Field code by R. Saavedra             !
! For more information please visit https://git.rsaavedra.xyz/?p=gcaf.git !
!---------------------------------------------------------------------!
!---------------------------------------------------------!
subroutine field(k)
!-----------------!
! External fields !
!-----------------!
!=========================================================!
	use shared,only:&
		keverg,& !keV to erg
		enaxo,&  !on axis en
		npert,&  !perturbation on/off
		nepot,&  !input pt profile pres/flat
		leno,&   !characteristic length
		ntmpi,&  !ions temp in keV
		ntmpe,&  !elec temp in keV
		pol,&    !poloidal flux
		polw,&   !on wall poloidal flux
		tht,&    !theta
		zet,&    !zeta
		rmin,&   !minor radius
		rmaj,&   !major radius
		cg,&     !poloidal current
		cgp,&    !g pol derivative
		ci,&     !toroidal current
		cip,&    !ri pol derivative
		qf,&     !q profile
		bf,&     !b field magnitude
		bm,&     !on axis b field magnitude
		dbdp,&   !b pol derivative
		dbdt,&   !b tht derivative
		dbdz,&   !b zet derivative
		pt,&     !electrostatic potential
		pto,&    !pt magnitude
		dptdp,&  !pt pol derivative
		dptdt,&  !pt tht derivative
		dptdz    !pt zet derivative
	implicit none
	integer(8):: k
	integer(8):: tmprat
	real(8):: pdum      !dummy pol
	real(8):: tdum      !dummy tht
	real(8):: zdum      !dummy zet
	real(8):: temp_prof !temp_prof function
	real(8):: tmpdum    !dummy temp
	real(8):: ew        !inverse aspect ratio
	real(8):: q0   		!on axis q profile
	real(8):: qw   		!on wall q profile
	real(8):: ra   		!r/a
	real(8):: rap  		!r/a pol derivative
	real(8):: lmsq 		!lambda squared
	real(8):: bfpft 	!b first term pol derivative
	real(8):: bfpst 	!b second term pol derivative
	real(8):: bt   		!b toroidal
	real(8):: bp   		!b poloidal
	real(8):: bpp  		!bp pol derivative
	real(8):: qfp  		!q pol derivative
	real(8):: dum
	real(8):: cipft
	real(8):: cipst
	real(8):: ciptt
!=========================================================!

	pdum= pol(k)
	tdum= tht(k)
	zdum= zet(k)
	ew= 1.d0*rmin/rmaj
	qw= 5.d0
	q0= 1.d0
	ra= sqrt(pdum/polw)
	rap= .5d0*ra/pdum
	lmsq= rmin*rmin*q0/(qw -q0)
	qf(k)= q0 +q0*ra*ra*rmin*rmin/lmsq
	qfp= q0*rmin*rmin/lmsq/polw
	cg(k)= rmaj/(1.d0 +ra*ew*cos(tdum))/leno/bm
	cgp(k)= rap*ew*cos(tdum)/(1.d0 +ra*ew*cos(tdum))
	ci(k)= rmin*ra*ew/(1.d0 +ra*ew*cos(tdum))/qf(k)/leno/bm
	cipft= rap*ew
	cipst= -qfp/qf(k)/qf(k)
	ciptt= rap*ew*cos(tdum)/(1.d0 +ra*ew*cos(tdum))
	cip(k)= cipft*cipst*ciptt

	!print *, cg(k)*leno*bm

	bt= bm
	bp= ra*ew*bm/qf(k)
	bpp= ew*bm*rap/qf(k) -ew*bm*ra*qfp/qf(k)/qf(k)
	bfpft= - ew*rap*cos(tdum)/(1.d0 +ew*ra*cos(tdum)) &
		/(1.d0 +ew*ra*cos(tdum))
	bfpst= bp*bpp/sqrt(bt*bt +bp*bp)
	bf(k)= sqrt(bt*bt +bp*bp)/bm/(1.d0 +ew*ra*cos(tdum))
	dbdp(k)= bfpst/bm/(1.d0 +ew*ra*cos(tdum)) &
		+bfpft*sqrt(bt*bt +bp*bp)/bm
	dbdt(k)= ew*ra*sin(tdum)*sqrt(bt*bt +bp*bp)/bm &
		/(1.d0 +ew*ra*cos(tdum))/(1.d0 +ew*ra*cos(tdum))
	dbdz(k)= 0.d0
	!simpler b
	!bf(k)= 1.d0 -ew*ra*cos(tdum)
	!dbdp(k)= -ew*rap*cos(tdum)
	!dbdt(k)= ew*ra*sin(tdum)
	!dbdz(k)= 0.d0

	pt(k)= 0.d0
	dptdp(k)= 0.d0
	dptdt(k)= 0.d0
	dptdz(k)= 0.d0
	dum= pdum/polw
	tmprat= ntmpi/ntmpi
	if(nepot.eq.0)then
		pt(k)= pto*dum*dum
		dptdp(k)= pto*dum/polw
		dptdt(k)= 0.d0
		dptdz(k)= 0.d0
	endif
	if(nepot.eq.1)then
		tmpdum= keverg/enaxo
		pt(k)= pto +tmpdum*temp_prof(dum,tmprat)
		dptdp(k)= pto -tmpdum*temp_prof(dum,tmprat)/polw
		dptdt(k)= 0.d0
		dptdz(k)= 0.d0
	endif
	if(npert.eq.1) call ptrba(k)
	if(pdum.ne.pdum) then
		write(0,'(a18)') 'POLOIDAL FLUX NaN'
		stop
	endif
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine ptrba(k)
!--------------------------------!
! Mynick alfa perturbation model !
! Phys. Fluids B 5, 1471 (1993)  !
!--------------------------------!
!=========================================================!
	use shared,only: &
		pi,&     !pi number
		namx,&   !input perturbation amplitude
		time,&   !time
		pol,&    !poloidal flux
		tht,&    !theta
		zet,&    !zeta
		al,&     !alfa perturbation
		dadp,&   !al pol derivative
		dadt,&   !al tht derivative
		dadz,&   !al zet derivative
		dadtim,& !al time derivative
		prtfrq,& !input perturbation frequency
		mmod,&   !input m mode number
		nmod,&   !input n mode number
		nphs,&	 !input perturbation phase
		polw     !on wall poloidal flux
	implicit none
	integer(8):: m   !m mode number
	integer(8):: n   !n mode number
	integer(8):: k   !particle index
	integer(8):: mi  !initial m mode
	integer(8):: mf  !final m mode
	real(8):: pdum   !pol dummy
	real(8):: tdum   !tht dummy
	real(8):: zdum   !zet dummy
	real(8):: ph     !phase
	real(8):: facx   !rmx/a
	real(8):: amx    !scaled namx
	real(8):: ra     !amp
	real(8):: rap    !amp prime
	real(8):: sn     !sine
	real(8):: cs     !cosine
	real(8):: rp     !Mynick model p=n for qnm>1
	real(8):: rotfrq !island rotation frequency
!=========================================================!
	al(k)= 0.d0
	dadp(k)= 0.d0
	dadt(k)= 0.d0
	dadz(k)= 0.d0
	dadtim(k)= 0.d0
	pdum= pol(k)
	tdum= tht(k)
	zdum= zet(k)

	ph= nphs*pi/180
	amx= namx
	mi=mmod
	mf=mmod
	n=nmod
	rp= 2 !n
	do m=mi,mf
		facx= 1.d0*m/(m +rp)
		ra= amx*(pdum/polw/facx)**m*(1.d0 -pdum/polw)**rp &
			/(1.d0 -facx)**rp
		rap= m*ra/pdum -rp*ra/polw/(1.d0 -pdum/polw)
		sn= sin(n*zdum -m*tdum +ph -prtfrq*time(k))
		cs= cos(n*zdum -m*tdum +ph -prtfrq*time(k))
		al(k)= al(k) +ra*sn
		dadp(k)= dadp(k) +rap*sn
		dadt(k)= dadt(k) -m*ra*cs
		dadz(k)= dadz(k) +n*ra*cs
		dadtim(k)= dadtim(k) -prtfrq*ra*cs
	enddo
	return
end subroutine
!---------------------------------------------------------!
