!---------------------------------------------------------------------!
! GCAF= Guiding Center Analytic Field code by R. Saavedra             !
! For more information please visit https://git.rsaavedra.xyz/?p=gcaf.git !
!---------------------------------------------------------------------!
!---------------------------------------------------------!
subroutine rcrd_log(n)
!-----------------!
! Record log file !
!-----------------!
!=========================================================!
	use shared
	implicit none
	integer(8):: n
	character(64):: arg
!=========================================================!

	if(n.eq.36) write(n,'(a8)') 'LOG FILE'
	write(n,'(a8)') ' '
	write(n,'(a10)') '# SWITCHES'
	if(ntype.eq.1) arg= 'traj'
	if(ntype.eq.2) arg= 'ense'
	if(ntype.eq.3) arg= 'poin'
	write(n,'(a26,a8)') 'simulation                ', adjustl(arg)
	write(n,'(a26,a4)') 'device                    ', 'tkmk'
	if(ncoll.eq.0) arg= 'off'
	if(ncoll.eq.1) arg= 'on'
	write(n,'(a26,a8)') 'collision                 ', adjustl(arg)
	if(npert.eq.0) arg= 'off'
	if(npert.eq.1) arg= 'on'
	write(n,'(a26,a8)') 'perturbation              ', adjustl(arg)
	write(n,'(a8)') ' '
	write(n,'(a12)') '# SIMULATION'
	write(arg,'(i8)') int(npart)
	write(n,'(a26,a8)') 'number_of_particles       ', adjustl(arg)
	write(arg,'(i8)') int(ntran)
	write(n,'(a26,a8)') 'time                      ', adjustl(arg)
	write(arg,'(i8)') int(nsdel)
	write(n,'(a26,a8)') 'step_delta                ', adjustl(arg)
	if(nfixs.eq.0) arg= 'off'
	if(nfixs.eq.1) arg= 'on'
	write(n,'(a26,a8)') 'step_fix                  ', adjustl(arg)
	write(arg,'(es10.2e2)') deler
	write(n,'(a26,a8)') 'max_error                 ', adjustl(arg)
	write(arg,'(i8)') int(ntraj)
	write(n,'(a26,a8)') 'trajectory_points         ', adjustl(arg)
	write(arg,'(i8)') int(ndisp)
	write(n,'(a26,a8)') 'dispersion_points         ', adjustl(arg)
	write(arg,'(i8)') int(ndist)
	write(n,'(a26,a8)') 'distribution_points       ', adjustl(arg)
	write(arg,'(es10.2e2)') r1
	write(n,'(a26,a8)') 'poincare_lower_bound      ', adjustl(arg)
	write(arg,'(es10.2e2)') r2
	write(n,'(a26,a8)') 'poincare_upper_bound      ', adjustl(arg)
	write(n,'(a8)') ' '
	write(n,'(a7)') '# FIELD'
	write(arg,'(i8)') int(bm)
	write(n,'(a26,a8)') 'field_magnitude           ', adjustl(arg)
	write(arg,'(i8)') int(rmaj)
	write(n,'(a26,a8)') 'major_radius              ', adjustl(arg)
	write(arg,'(i8)') int(rmin)
	write(n,'(a26,a8)') 'minor_radius              ', adjustl(arg)
	write(arg,'(es10.2e2)') namx
	write(n,'(a26,a8)') 'perturbation_amplitude    ', adjustl(arg)
	write(arg,'(i8)') int(nfrq)
	write(n,'(a26,a8)') 'perturbation_frequency    ', adjustl(arg)
	write(arg,'(i8)') int(nphs)
	write(n,'(a26,a8)') 'perturbation_phase        ', adjustl(arg)
	write(arg,'(i8)') int(nmod)
	write(n,'(a26,a8)') 'n_mode                    ', adjustl(arg)
	write(arg,'(i8)') int(mmod)
	write(n,'(a26,a8)') 'm_mode                    ', adjustl(arg)
	write(arg,'(es10.2e2)') npto
	write(n,'(a26,a8)') 'electrostatic_potential   ', adjustl(arg)
	if(nepot.eq.0) arg= 'cons'
	if(nepot.eq.1) arg= 'pres'
	write(n,'(a26,a8)') 'potential_profile         ', adjustl(arg)
	write(n,'(a8)') ' '
	write(n,'(a8)') '# PLASMA'
	if(nspec.eq.0) arg= 'ions'
	if(nspec.eq.1) arg= 'elec'
	if(nspec.eq.2) arg= 'both'
	if(nspec.eq.0) arg= 'ptch'
	if(nspec.eq.1) arg= 'ener'
	if(nspec.eq.2) arg= 'both'
	!write(n,'(a26,a8)') 'scattering_operator       ', adjustl(arg)
	write(arg,'(i8)') mi
	write(n,'(a26,a8)') 'ions_mass                 ', adjustl(arg)
	write(arg,'(i8)') zi
	write(n,'(a26,a8)') 'ions_charge               ', adjustl(arg)
	write(arg,'(i8)') ntmpi
	write(n,'(a26,a8)') 'iones_temperature         ', adjustl(arg)
	write(arg,'(i8)') ntmpe
	write(n,'(a26,a8)') 'electrons_temperature     ', adjustl(arg)
	if(ntemp.eq.0) arg= 'flat'
	if(ntemp.eq.1) arg= 'gaus'
	write(n,'(a26,a8)') 'temperature_profile       ', adjustl(arg)
	write(arg,'(es10.2e2)') deno
	write(n,'(a26,a8)') 'density                   ', adjustl(arg)
	if(ndens.eq.0) arg= 'flat'
	if(ndens.eq.1) arg= 'para'
	write(n,'(a26,a8)') 'density_profile           ', adjustl(arg)
	write(n,'(a8)') ' '
	write(n,'(a15)') '# TEST PARTICLE'
	if(cs.eq.1) arg= 'ion'
	if(cs.eq.-1) arg= 'ele'
	write(n,'(a26,a8)') 'particle_type             ', adjustl(arg)
	write(arg,'(i8)') mt
	write(n,'(a26,a8)') 'mass                      ', adjustl(arg)
	write(arg,'(i8)') zt
	write(n,'(a26,a8)') 'charge                    ', adjustl(arg)
	write(arg,'(i8)') ptchgrad
	write(n,'(a26,a8)') 'pitch_angle               ', adjustl(arg)
	write(arg,'(es10.2e2)') nengo
	write(n,'(a26,a8)') 'energy                    ', adjustl(arg)
	if(nends.eq.0) arg= 'mon'
	if(nends.eq.1) arg= 'max'
	write(n,'(a26,a8)') 'energy_distribution       ', adjustl(arg)
	write(arg,'(i8)') nrado
	write(n,'(a26,a8)') 'radius                    ', adjustl(arg)
	write(arg,'(i8)') thtgrad
	write(n,'(a26,a8)') 'theta                     ', adjustl(arg)
	write(arg,'(i8)') zetgrad
	write(n,'(a26,a8)') 'zeta                      ', adjustl(arg)
	write(n,'(a8)') ' '
	write(n,'(a14)') '# OUTPUT FILES'
	write(arg,'(a64)') trim(lngargtraj)
	write(n,'(a64)') adjustl(arg)
	write(arg,'(a64)') trim(lngargpoin)
	write(n,'(a64)') adjustl(arg)
	write(arg,'(a64)') trim(lngargdist)
	write(n,'(a64)') adjustl(arg)
	write(arg,'(a64)') trim(lngargdisp)
	write(n,'(a64)') adjustl(arg)
	write(arg,'(a64)') trim(lngargcoef)
	write(n,'(a64)') adjustl(arg)
	write(n,'(a8)') ' '
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine rcrd_traj(k,nskip)
!-------------------!
! Record trajectory !
!-------------------!
!=========================================================!
	use shared,only: &
		pi,&     !pi number
		keverg,& !keV to erg
		rc,&     !speed of light
		ptch,&   !pitch
		en,&     !energy
		pol,&    !poloidal flux
		tht,&    !theta
		zet,&    !zeta
		rho,&    !parallel gyroradius
		time,&   !time variable
		polw,&   !on wall poloidal flux
		leno,&   !characteristic length
		omego,&  !characteristic time freq
		bf,&     !b field
		mt,&     !test particle mass in mp
		enaxo,&  !on axis en
		pn01,&
		pn02,&
		pn03,&
		pn04,&
		pn05,&
		pn06,&
		pn07,&
		pn08,&
		pn09,&
		pn10,&
		trgttraj,& !traj targets
		rcrdtraj !record traj points
	implicit none
	integer(8):: k
	integer(8):: i
	integer(8):: mdum
	integer(8):: npts
	integer(8):: nskip
	real(8):: tdum  !dummy variables
	real(8):: xdum
	real(8):: zdum
	real(8):: xxdum
	real(8):: yydum
	real(8):: zzdum
	real(8):: xxproj !proj functions
	real(8):: yyproj
	real(8):: zzproj
	real(8):: bigr   !cylindrical coordinates
	real(8):: bigz
	real(8):: bigp
	real(8):: pparl  !parallel and perpendicular moments
	real(8):: pperp
!=========================================================!

	!tokamak stellarator record
	tdum= time(k)/omego
	xdum= xxproj(pol(k),tht(k),0.d0)
	zdum= zzproj(pol(k),tht(k),0.d0)
	xxdum= xxproj(pol(k),tht(k),zet(k))
	yydum= yyproj(pol(k),tht(k),zet(k))
	zzdum= zzproj(pol(k),tht(k),zet(k))
	bigr= sqrt(xxdum*xxdum +yydum*yydum)
	bigz= zzdum

	!bound 0< zet< 2*pi
	mdum= zet(k)/2/pi
	bigp= zet(k) -2*pi*mdum
	if(bigp.lt.0.d0) bigp= bigp +2*pi

	pparl= mt*leno*omego*rho(k)*bf(k)/rc
	pperp= mt*leno*omego*sqrt(en(k) -.5d0*rho(k)*rho(k)*bf(k)*bf(k))/rc
	if(pperp.ne.pperp) pperp= 0.d0

	pn01= tdum
	pn02= 2*sqrt(pol(k))*leno
	pn03= tht(k)
	pn04= zet(k)
	pn05= pparl
	pn06= pperp
	pn07= ptch(k)
	pn08= en(k)*enaxo/keverg
	pn09= bigr
	pn10= zzdum

	do i=1,size(trgttraj)
		rcrdtraj(nskip,i)= trgttraj(i)
	enddo
	nskip= nskip +1
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine rcrd_coef(ind,slp1,slp2,slp3)
!------------------------------------------------------!
! Obtain diffusion coefficients from linear regression !
! of the previously recorded values for the dispersion !
!------------------------------------------------------!
!=========================================================!
	use shared,only: &
		npart,&  !input number of particles
		ndisp,&  !input number of disp points
		ndist,&  !input number of dist points
		omego,&  !characteristic time frequency
		dft,&    !time for disp
		rsqav,&  !rad sq av
		rav,&    !rad av
		ensqav,& !en sq av
		enav,&   !en av
		pzsqav,& !pz sq av
		pzav,&   !pz av
		fm       !output format
	implicit none
	integer(8):: &
		j,k,&  !indices
		trun,& !truncate time
		ind    !snapshot index
	real(8):: &
		tdum,& !time dummy
		dff1,& !r dispersion
		dff2,& !en dispersion
		dff3,& !pz dispersion
		av,&   !average
		sqav,& !squared average
		slp1,& !r slope
		slp2,& !en slope
		slp3,& !pz slope
		slpdum,& !slope dummy
		devdum
	real(8),dimension(ndist,ndisp):: &
		tmat,&
		dff1mat,&
		dff2mat,&
		dff3mat
	real(8),dimension(:),allocatable:: &
		tmatdum,&
		dff1matdum,&
		dff2matdum,&
		dff3matdum
!=========================================================!
	trun= ind*ndisp/ndist
	do j=1,trun
		tdum= dft(j)/omego
		tmat(ind,j)= tdum
		sqav= rsqav(j)/npart
		av= rav(j)/npart
		dff1= sqav- av*av
		dff1mat(ind,j)= dff1
		sqav= ensqav(j)/npart
		av= enav(j)/npart
		dff2= sqav- av*av
		dff2mat(ind,j)= dff2
		sqav= pzsqav(j)/npart
		av= pzav(j)/npart
		dff3= sqav- av*av
		dff3mat(ind,j)= dff3
	enddo
	allocate(&
		tmatdum(trun),&
		dff1matdum(trun),&
		dff2matdum(trun),&
		dff3matdum(trun))
	tmatdum= tmat(ind,1:trun)
	dff1matdum= dff1mat(ind,1:trun)
	dff2matdum= dff2mat(ind,1:trun)
	dff3matdum= dff3mat(ind,1:trun)
	call least_square(trun,tmatdum,dff1matdum,slpdum,devdum)
	slp1= slpdum/2
	call least_square(trun,tmatdum,dff2matdum,slpdum,devdum)
	slp2= slpdum/2
	call least_square(trun,tmatdum,dff3matdum,slpdum,devdum)
	slp3= slpdum/2
	deallocate(&
		tmatdum,&
		dff1matdum,&
		dff2matdum,&
		dff3matdum)
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine write_traj
!------------------!
! Write trajectory !
!------------------!
!=========================================================!
	use shared,only: &
		fm,&       !output format
		ntraj,&    !input number of traj points
		trgttraj,& !traj targets
		rcrdtraj   !record trajectory points
	implicit none
	integer(8):: i
!=========================================================!

	do i=1,ntraj
		write(43,fm) rcrdtraj(i,1:size(trgttraj))
	enddo
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine write_dist
!----------------------!
! Record distributions !
!----------------------!
!=========================================================!
	use shared,only: &
		npart,&      !input number of particles
		ndist,&      !input number of dist points
		ndisp,&      !input number of disp points
		distt,&      !time distribution
		distr,&      !radial distribution
		diste,&      !en distribution
		distp,&      !pitch angle distribution
		trgtdist,&   !disp targets,&
		lngargdist,& !dist.plt user string
		pn41,&
		pn42,&
		pn43,&
		fm           !output format
	implicit none
	integer(8):: j
	integer(8):: k
	integer(8):: n
	real(8):: coef1
	real(8):: coef2
	real(8):: coef3
	character(32):: label
	character(32):: arg
!=========================================================!
	do j=1,ndist
		n= 60+j
		call rcrd_coef(j,coef1,coef2,coef3)
		write(arg,'(i2)') j
		write(n,*) '# dist_'//adjustl(trim(arg))//&
			'    time rcoef encoef pzcoef'
		write(n,fm) distt(j),coef1,coef2,coef3
		write(n,*) '# ',trim(lngargdist)
		do k=1,npart
			pn41= distr(k,j)
			pn42= diste(k,j)
			pn43= distp(k,j)
			write(n,fm) trgtdist(:)
		enddo
	enddo
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine write_disp
!---------------------------!
! Write statistical moments !
!---------------------------!
!=========================================================!
	use shared,only: &
		npart,&  !input number of particles
		ndisp,&  !input number of disp points
		npto,&   !input electrostatic potential in keV
		omego,&  !characteristic time freq
		dft,&    !time for disp
		rav,&    !rad av
		rsqav,&  !rad sq av
		rcuav,&  !rad cu av
		rfoav,&  !rad fo av
		enav,&   !en av
		ensqav,& !en sq av
		encuav,& !en cu av
		enfoav,& !en fo av
		pzav,&   !pz av
		pzsqav,& !pz sq av
		pzcuav,& !pz cu av
		pzfoav,& !pz fo av
		deno,&   !plasma density in cm^-3
		pol,&    !poloidal flux
		tht,&    !theta
		pn21,&   !tim disp
		pn22,&   !dsr disp
		pn23,&   !dse disp
		pn24,&   !dsz disp
		pn25,&   !skr disp
		pn26,&   !ske disp
		pn27,&   !skz disp
		pn28,&   !flr disp
		pn29,&   !fle disp
		pn30,&   !flz disp
		pn31,&   !den coef
		pn32,&   !pot coef
		pn33,&   !flx coef
		pn34,&   !dfr coef
		pn35,&   !dfe coef
		pn36,&   !dfp coef
		pn37,&   !sdr coef
		pn38,&   !sde coef
		pn39,&   !sdp coef
		trgtdisp,& !disp targets
		trgtcoef,& !coef targets
		flux,&     !particle flux
		fm         !output format
	implicit none
	integer(8):: j,k
	real(8):: tdum !time dummy
	real(8):: xdum !x dummy
	real(8):: zdum !z dummy
	real(8):: var1 !r dispersion
	real(8):: var2 !en dispersion
	real(8):: var3 !pz dispersion
	real(8):: skw1 !r skewness
	real(8):: skw2 !en skewness
	real(8):: skw3 !pz skewness
	real(8):: fla1 !r flatness
	real(8):: fla2 !en flatness
	real(8):: fla3 !pz flatness
	real(8):: sdev !standard deviation
	real(8):: av   !average
	real(8):: sqav !square average
	real(8):: cuav !cubic average
	real(8):: foav !fourth average
	real(8):: flx1 !flux/npart
	real(8):: slp1 !r slope
	real(8):: slp2 !en slope
	real(8):: slp3 !pz slope
	real(8):: dev1 !standard deviation 1
	real(8):: dev2 !standard deviation 2
	real(8):: dev3 !standard deviation 3
	real(8):: xxproj !x proj fun
	real(8):: zzproj !y proj fun
	real(8),dimension(ndisp):: tvec    !time vector
	real(8),dimension(ndisp):: var1vec !variance r vector
	real(8),dimension(ndisp):: var2vec !variance e vector
	real(8),dimension(ndisp):: var3vec !variance pz vector
!=========================================================!
	do j=1,ndisp
		tdum= dft(j)/omego
		av= rav(j)/npart
		sqav= rsqav(j)/npart
		cuav= rcuav(j)/npart
		foav= rfoav(j)/npart
		sdev= sqrt(abs(sqav- av*av))
		var1= sqav- av*av
		skw1= (cuav -3*av*sqav +2*av*av*av)/sdev/sdev/sdev
		fla1= (foav -4*av*cuav +6*sqav*av*av -3*av*av*av*av)&
			/sdev/sdev/sdev/sdev
		av= enav(j)/npart
		sqav= ensqav(j)/npart
		cuav= encuav(j)/npart
		foav= enfoav(j)/npart
		sdev= sqrt(abs(sqav- av*av))
		var2= sqav- av*av
		skw2= (cuav -3*av*sqav +2*av*av*av)/sdev/sdev/sdev
		fla2= (foav -4*av*cuav +6*sqav*av*av -3*av*av*av*av)&
			/sdev/sdev/sdev/sdev
		av= pzav(j)/npart
		sqav= pzsqav(j)/npart
		cuav= pzcuav(j)/npart
		foav= pzfoav(j)/npart
		sdev= sqrt(abs(sqav- av*av))
		var3= sqav- av*av
		skw3= (cuav -3*av*sqav +2*av*av*av)/sdev/sdev/sdev
		fla3= (foav -4*av*cuav +6*sqav*av*av -3*av*av*av*av)&
			/sdev/sdev/sdev/sdev
		pn21= tdum
		pn22= var1
		pn23= var2
		pn24= var3
		pn25= skw1
		pn26= skw2
		pn27= skw3
		pn28= fla1
		pn29= fla2
		pn30= fla3
		write(50,fm) trgtdisp(:)
		tvec(j)= tdum
		var1vec(j)= var1
		var2vec(j)= var2
		var3vec(j)= var3
	enddo

	call least_square(ndisp,tvec,var1vec,slp1,dev1)
	call least_square(ndisp,tvec,var2vec,slp2,dev2)
	call least_square(ndisp,tvec,var3vec,slp3,dev3)
	flx1= 1.d0*flux/npart

	pn31= deno
	pn32= npto
	pn33= flx1
	pn34= slp1
	pn35= slp2
	pn36= slp3
	pn37= dev1
	pn38= dev2
	pn39= dev3
	write(70,'(i8,12es12.4e2)') npart,trgtcoef(:)
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine write_poin(k)
!-----------------------------------!
! Record kinetic Poincaré plot data !
!-----------------------------------!
!=========================================================!
	use shared,only: &
		pi,&       !pi number
		time,&     !time
		pol,&      !poloidal flux
		tht,&      !theta
		zet,&      !zeta
		rho,&      !parallel gyroradius
		polw,&     !on wall poloidal flux
		prtfrq,&   !input perturbation frequency
		zplo,&     !poin zeta plane
		nmod,&     !input n mode number
		rmin,&     !minor radius
		leno,&     !charasteristic length
		r1,&       !poin lower bound
		r2,&       !poin upper bound
		pold,&     !old pol for poin
		told,&     !old tht for poin
		zold,&     !old zet for poin
		pn11,& 	   !rad poin
		pn12,&     !tht poin
		pn13,&     !rcy poin
		pn14,&     !zcy poin
		trgtpoin,& !poin targets
		fm         !output format
	implicit none
	integer(8):: k    !index
	integer(8):: mdum !dummy bound
	real(8):: pdum    !pol dummy
	real(8):: rdum    !rad dummy
	real(8):: tdum    !tht dummy
	real(8):: zdum    !zet dummy
	real(8):: dum     !zold dummy
	real(8):: xxproj  !x proj fun
	real(8):: yyproj  !y proj fun
	real(8):: zzproj  !z proj fun
	real(8):: xpp     !x proj
	real(8):: ypp     !y proj
	real(8):: zpp     !z proj
	real(8):: rpp     !r cyl
	real(8):: zpl     !z cyl
	real(8):: del     !delta
	real(8):: p1      !pol lower bound
	real(8):: p2      !pol upper bound
	real(8):: r1dum   !pol lower bound dummy
	real(8):: r2dum   !pol upper bound dummy
!=========================================================!
	r1dum= r1*rmin
	r2dum= r2*rmin
	p1= .5d0*r1dum*r1dum/leno/leno
	p2= .5d0*r2dum*r2dum/leno/leno
	del= pi/2000
	zpl= zplo
	if(zpl.eq.0.d0) zpl= pi/100

	!bound 0< zet< 2*pi
	mdum= (nmod*zet(k) -prtfrq*time(k))/2/pi
	zdum= (nmod*zet(k) -prtfrq*time(k)) -2*pi*mdum
	if(zdum.lt.0.d0) zdum= zdum +2*pi

	!check if particle has crossed plane
	dum= zold(k)
    if(zpl.gt.0.d0) then
    	if(dum.gt.zpl+del) go to 311
    	if(dum.lt.zpl) go to 311
    endif
    if(zpl.lt.0.d0) then
        if(dum.lt.-zpl-del) go to 311
        if(dum.gt.-zpl) go to 311
    endif

	!linear interpolate and evaluate at plane
	pdum= (pold(k)*(zdum -zpl) +pol(k)*(zpl -zold(k)))/(zdum -zold(k))
	tdum= (told(k)*(zdum -zpl) +tht(k)*(zpl -zold(k)))/(zdum -zold(k))

	!bound p1< pol< p2
	if(pdum.lt.p1) go to 311
	if(pdum.gt.p2) go to 311
	xpp= xxproj(pdum,tdum,zpl)
	ypp= yyproj(pdum,tdum,zpl)
	zpp= zzproj(pdum,tdum,zpl)
	rpp= sqrt(xpp*xpp +ypp*ypp)

	!bound -pi< thet< pi
	mdum= tdum/2/pi
	tdum= tdum -2*pi*mdum
	if(tdum.lt.-pi) tdum= tdum +2*pi
	if(tdum.gt.pi) tdum= tdum -2*pi
	rdum= sqrt(pdum/polw)
	pdum= pdum/polw

	pn11= rdum
	pn12= tdum
	pn13= rpp
	pn14= zpp
	write(48,'(i8,8es12.4e2)') k,trgtpoin(:)

	311 continue
	pold(k)= pol(k)
	told(k)= tht(k)
	zold(k)= zdum
	return
end subroutine
!---------------------------------------------------------!
