!---------------------------------------------------------------------!
! GCAF= Guiding Center Analytic Field code by R. Saavedra             !
! Email of the author: saavestro@gmail.com                            !
! For more information please visit https://github.com/Saavestro/gcaf !
!---------------------------------------------------------------------!
!---------------------------------------------------------!
subroutine scatr(k)
!-----------------------------------------------------!
! Boozer and Kuo-Petravic Phys. Fluids 24, 851 (1981) !
!-----------------------------------------------------!
!=========================================================!
	use shared,only: &
		pi,&     !pi number
		keverg,& !keV to erg
		rq,&     !elec charge
		mp,&     !proton mass
		ndens,&  !dens profile flat/para
		ntemp,&  !temp profile flat/gaus
		ntmpi,&	 !input ions temp in keV
		ntmpe,&	 !input elec temp in keV
		deno,&   !plasma dens in cm^-3
		ptch,&   !pitch angle
		en,&     !energy
		pol,&    !poloidal flux
		rho,&    !parallel gyroradius
		polw,&   !on wall poloidal flux
		zi,&	 !ions charge in protons charge
		zt,&     !test charge in protons charge
		mt,&     !test particle mass in protons mass
		mi,&     !ions mass in protons mass
		enaxo,&  !on axis particle en 
		leno,&   !characteristic length
		omego,&  !characteristic time freq
		dt,&     !time step
		bf,&     !b field 
		pt,&     !electrostatic potential
		mu       !b moment
	implicit none
	integer(4):: k
	real(8):: dni       !ions density in cm-3
	real(8):: ti        !ions temperature in code units
	real(8):: te        !electrons temperature code units
	real(8):: vti       !ions thermal velocity
	real(8):: dum
	real(8):: fi
	real(8):: gi
	real(8):: efaci
	real(8):: eface
	real(8):: psig      !pitch operator random sign
	real(8):: esig      !en operator random sign
	real(8):: xd		!pol/polw for dens_prof
	real(8):: xt		!pol/polw for temp_prof
	real(8):: vprt      !test particle thermal velocity
	real(8):: dens_prof !dens profile function
	real(8):: temp_prof !temp profile function
	real(8):: ranx      !random function
	real(8):: vte       !electron thermal velocity
	real(8):: fe
	real(8):: ge
	real(8):: dumi      !vprt/vti
	real(8):: dume      !vprt/vte
	real(8):: clogi     !ions Coulomb log
	real(8):: cloge     !elec Coulomb log
	real(8):: tikev     !ions temp in KeV
	real(8):: tekev     !elec temp in KeV
	real(8):: dne       !elec dens in cm-3
	real(8):: eekev     !test kinetic en in keV
	real(8):: ekin      !test kinetic en
	real(8):: ekmin     !test en minimum
	real(8):: colbi     !ions braginskii freq
	real(8):: colbe     !elec braginskii freq
	real(8):: coldi     !ions deflection freq
	real(8):: colde     !elec deflection freq
	real(8):: colei     !ions en freq
	real(8):: colee     !elec en freq
	real(8):: colds     !deflection freq sum
	real(8):: coles     !en freq sum
	real(8):: velo      !characteristic velocity
	real(8):: clogtest 
	real(8):: bmintest
	real(8):: bmaxtest
	real(8):: mutest
	real(8):: usqtest
!=========================================================!

	!SET PARAMETERS
	xd= pol(k)/polw
	if(ndens.eq.0) xd= 0.d0
	xt= pol(k)/polw
	if(ntemp.eq.0) xt= 0.d0
	ekin= en(k)
	eekev= enaxo*en(k)/keverg 
	dni= dens_prof(xd,deno)
	dne= dni*zi 
	tikev= temp_prof(xt,ntmpi/ntmpi)
	tekev= temp_prof(xt,ntmpe/ntmpi)
	ti= keverg*tikev/enaxo 
	te= keverg*tekev/enaxo 
	ekmin= 1.d-3*ti
	vprt= sqrt(2*ekin/mt)
	vti= sqrt(ti/mi)
	vte= sqrt(1836*te) 
	velo= leno*omego
	
	!IONS
	dumi= vprt/vti
	fi= ((1.d0- .5d0/dumi/dumi)*erf(dumi) &
		+exp(-dumi*dumi)/sqrt(pi)/dumi)/dumi/dumi/dumi
	gi = .5d0*(erf(dumi) -2.d0*dumi*exp(-dumi**2)/sqrt(pi))/dumi**3 !en
	efaci= ekin -exp(-dumi**2)*ti/gi/sqrt(pi) !en
	clogi= 23.d0 -log(zt*zt*zt*sqrt(2*dni)/(1.d+3*tikev)**1.5d0)
	!!
	mutest= mp*mi*mt/(mi +mt)
	usqtest= leno*leno*omego*omego*vprt*vprt
	bmintest= zt*zi*rq*rq/mutest/usqtest
	bmaxtest= 1/sqrt(dni*zi*zi*rq*rq/tikev/keverg +dne*rq*rq/tekev/keverg)
	clogtest= log(bmaxtest/bmintest)
	!!
	colbi= 1.2717d+11*(1.d0*zi*zi*zi/mi/mi)*dni*clogi &
		/velo/velo/velo/vti/vti/vti
	coldi= abs(1.8799d0*(zt*zt*mi*mi/zi/zi*mt*mt)*colbi*fi)/omego
	colei= abs(3.7599d0*(zt**2*mi/zi**2/mt)*colbi*gi)/omego !en

	!ELECTRONS
	dume= vprt/vte
	fe= ((1.d0- .5d0/dume**2)*erf(dume) &
		+exp(-dume**2)/sqrt(pi)/dume)/dume**3
	ge= .5d0*(erf(dume) -2.d0*dume*exp(-dume**2)/sqrt(pi))/dume**3 !en
	eface= ekin -exp(-dume**2)*te/ge/sqrt(pi) !en
	if(1.d+3*tekev.le. 5.46d-4*zt*1.d+3*tikev/mt) &
		cloge= 23.d0 -log(zt*sqrt(dne)/(1.d+3*tekev)**1.5d0)
	if(1.d+3*tekev.gt. 5.46d-4*zt*1.d+3*tikev/mt .and. &
		1.d+3*tekev.le. 10.d0*zt**2) &
		cloge= 30.d0 -log(zt*sqrt(dne)/(1.d+3*tekev)**1.5d0/mt)
	if(1.d+3*tekev.gt. 10.d0*zt**2) &
		cloge= 24.d0 -log(sqrt(dne)/1.d+3/tekev)
	colbe= 6.0634d+17*zi**2*dne*cloge/(leno*omego)**3/vte**3
	colde= abs(1.3293d0*(zt**2/zi/mi**2)*colbe*fe/1836**2)/omego
	colee= abs(2.6586d0*(zt**2/zi**2/mt)*colbe*ge/1836)/omego !en

	!FREQUENCIES
	colds= coldi +colde
	coles= ti*colei +te*colee !en
	
	!COLLISION STEP
	!pitch-angle step 
	dum= max(0.d0,1.d0- ptch(k)**2)
	psig= sign(1.d0,ranx()-.5d0)*sqrt(dum*colds*dt(k))
	ptch(k)= ptch(k)*(1.d0 -colds*dt(k))
	ptch(k)= ptch(k) +psig
	ptch(k)= min(1.d0,ptch(k))
	ptch(k)= max(-1.d0,ptch(k))
	!en step
	esig= sign(1.d0,ranx()-.5d0)*2.d0*sqrt(coles*ekin*dt(k))
	ekin= ekin -2.d0*colei*efaci*dt(k) -2.d0*colee*eface*dt(k)
	ekin= ekin +esig
	ekin= max(ekin,ekmin)
	!en(k)= ekin 
	
	!change rho, mu, due to scattering
	rho(k)= ptch(k)*sqrt(2*en(k))/bf(k)
	mu(k)= en(k)/bf(k)- .5d0*rho(k)*rho(k)*bf(k)
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
function dens_prof(xd,den)
!=========================================================!
	implicit none
	real(8):: dens_prof !dens_prof function
	real(8):: den       !input density
	real(8):: xd        !variable
!=========================================================!
	dens_prof= den*(1.d0- xd*xd)
	return
end function
!---------------------------------------------------------!


!---------------------------------------------------------!
function temp_prof(xd,temp)
!=========================================================!
	use shared,only: &
		pi
	implicit none
	integer(8):: temp   !input temperature	
	real(8):: temp_prof !temp_prof function
	real(8):: xd     	!variable
	real(8):: kern    	!gaussian kernel
!=========================================================!
	kern= xd
	temp_prof= temp*exp(-kern)
	return
end function
!---------------------------------------------------------!
