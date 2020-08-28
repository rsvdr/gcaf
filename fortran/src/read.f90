!---------------------------------------------------------------------!
! GCAF= Guiding Center Analytic Field code by R. Saavedra             !
! Email of the author: saavestro@gmail.com                            !
! For more information please visit https://github.com/Saavestro/gcaf !
!---------------------------------------------------------------------!
!---------------------------------------------------------!
subroutine cmdlin
!------------------------!
! command line arguments !
!------------------------!
!=========================================================!
	use shared
	implicit none
	integer(4):: i         !index
	real(8):: narg         !integer argumet
	character(64):: arg    !real argument
	character(64):: dir    !directory
	character(64):: outdir !output directory
!=========================================================!

	i= 0
	do
		call get_command_argument(i,arg)
		if(len_trim(arg).eq.0) exit
		
		if(trim(arg).eq.'-i')then 
			call get_command_argument(i+1,arg)
			open(32,file=adjustl(trim(arg)),status='old',action='read')
			call read_data
		endif
		if(trim(arg).eq.'-o')then 
			call get_command_argument(i+1,outdir)
		endif		
		if(trim(arg).eq.'-sim') then
			call get_command_argument(i+1,arg)
			if(adjustl(trim(arg)).eq.'traj') ntype= 1
			if(adjustl(trim(arg)).eq.'ense') ntype= 2
			if(adjustl(trim(arg)).eq.'poin') ntype= 3
		endif
		if(trim(arg).eq.'-col') then
			call get_command_argument(i+1,arg)
			if(adjustl(trim(arg)).eq.'on') ncoll= 1
			if(adjustl(trim(arg)).eq.'off') ncoll= 0
		endif
		if(trim(arg).eq.'-prt') then
			call get_command_argument(i+1,arg)
			if(adjustl(trim(arg)).eq.'on') npert= 1
			if(adjustl(trim(arg)).eq.'off') npert= 0
		endif
		if(trim(arg).eq.'-scr') then 
			call get_command_argument(i+1,arg)
			if(adjustl(trim(arg)).eq.'on') nscrn= 1
			if(adjustl(trim(arg)).eq.'off') nscrn= 0
		endif
		if(trim(arg).eq.'-npr') then
			call get_command_argument(i+1,arg)
			read(arg,*) narg
			npart= narg
		endif
		if(trim(arg).eq.'-pch') then
			call get_command_argument(i+1,arg)
			read(arg,*) narg
			ptchgrad= narg
		endif
		if(trim(arg).eq.'-eng') then
			call get_command_argument(i+1,arg)
			read(arg,*) narg
			nengo= narg
		endif
		if(trim(arg).eq.'-amp') then
			call get_command_argument(i+1,arg)
			read(arg,*) narg
			namx= narg
		endif
		if(trim(arg).eq.'-pot') then 
			call get_command_argument(i+1,arg)
			read(arg,*) narg
			npto= narg
		endif
		if(trim(arg).eq.'-den') then
			call get_command_argument(i+1,arg)
			read(arg,*) narg
			deno= narg
		endif
		if(trim(arg).eq.'-tht') then
			call get_command_argument(i+1,arg)
			read(arg,*) narg
			thtgrad= narg
		endif
		i= i +1
	enddo	
	call open_files(outdir)
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine read_data
!----------------------------!
! Reads data from input file !
!----------------------------!
!=========================================================!
	use shared
	implicit none
	integer(8):: i       !index
	integer(8):: eof     !equation of 
	character(32):: text !text line
	character(32):: arg  !argument
	character(64):: line !line
	character(64):: strp !strip function
!=========================================================!
	
	ntype= 1
	nconf= 0
	ncoll= 0
	npert= 0
	npart= 1
	ntran= 80
	nsdel= 1000
	nfixs= 0
	deler= 1.d-6
	ntraj= 10000
	ndisp= 400
	ndist= 4
	r1= 0.2d0
	r2= 0.8d0
	bm= 1
	rmaj= 150
	rmin= 50
	namx= 1.d0
	nfrq= 0
	mmod= 2
	nmod= 1
	nphs= 0
	pto= 0
	nepot= 0
	nspec= 2
	nscat= 0
	cs= 1
	mi= 1
	zi= 1
	ntmpi= 1
	ntmpe= 1
	deno= 2.d+14
	mt= 1
	zt= 1
	ptchgrad= 0
	nengo= 1
	nends= 0
	nrado= rmin/2
	thtgrad= 0
	zetgrad= 0
	lngargtraj= 'traj.plt    tim rad tht zet par per eng pch rcy zcy'
	lngargpoin= 'poin.plt    rad tht rcy zcy'
	lngargdist= 'dist.plt    rad eng pch'
	lngargdisp= 'disp.plt    tim pol eng pzt'
	lngargcoef= 'coef.plt    den dfr dfe dfp sdr sde sdp'
	do
		read(32,'(a)',iostat=eof) line
		text=''
		if(line(1:1).ne.'#'.and.len_trim(line).ne.0) read(line,*) text,arg
		if(trim(text).eq.'simulation') then 
			if(adjustl(trim(arg)).eq.'traj') ntype= 1
			if(adjustl(trim(arg)).eq.'ense') ntype= 2
			if(adjustl(trim(arg)).eq.'poin') ntype= 3
		endif
		if(trim(text).eq.'device') then 
			if(adjustl(trim(arg)).eq.'tkmk') nconf= 0
			if(adjustl(trim(arg)).eq.'stel') nconf= 1
		endif
		if(trim(text).eq.'collisions') then 
			if(adjustl(trim(arg)).eq.'on') ncoll= 1
			if(adjustl(trim(arg)).eq.'off') ncoll= 0
		endif
		if(trim(text).eq.'perturbation') then 
			if(adjustl(trim(arg)).eq.'on') npert= 1
			if(adjustl(trim(arg)).eq.'off') npert= 0
		endif
		if(trim(text).eq.'number_of_particles') read(arg,*) npart
		if(trim(text).eq.'time') read(arg,*) ntran
		if(trim(text).eq.'step_delta') read(arg,*) nsdel
		if(trim(text).eq.'step_fix') then 
			if(adjustl(trim(arg)).eq.'on') nfixs= 1
			if(adjustl(trim(arg)).eq.'off') nfixs= 0
		endif
		if(trim(text).eq.'max_error') read(arg,*) deler
		if(trim(text).eq.'trajectory_points') read(arg,*) ntraj
		if(trim(text).eq.'dispersion_points') read(arg,*) ndisp
		if(trim(text).eq.'distribution_points') read(arg,*) ndist
		if(trim(text).eq.'poincare_lower_bound') read(arg,*) r1
		if(trim(text).eq.'poincare_upper_bound') read(arg,*) r2
		if(trim(text).eq.'field_magnitude') read(arg,*) bm
		if(trim(text).eq.'major_radius') read(arg,*) rmaj
		if(trim(text).eq.'minor_radius') read(arg,*) rmin
		if(trim(text).eq.'perturbation_amplitude') read(arg,*) namx
		if(trim(text).eq.'perturbation_frequency') read(arg,*) nfrq
		if(trim(text).eq.'perturbation_phase') read(arg,*) nphs
		if(trim(text).eq.'m_mode') read(arg,*) mmod
		if(trim(text).eq.'n_mode') read(arg,*) nmod
		if(trim(text).eq.'electrostatic_potential') read(arg,*) npto
		if(trim(text).eq.'potential_profile')then
			if(adjustl(trim(arg)).eq.'flat') nepot= 0
			if(adjustl(trim(arg)).eq.'pres') nepot= 1
		endif
		if(trim(text).eq.'particle_species')then
			if(adjustl(trim(arg)).eq.'ions') nspec= 0
			if(adjustl(trim(arg)).eq.'elec') nspec= 1
			if(adjustl(trim(arg)).eq.'both') nspec= 2
		endif
		if(trim(text).eq.'scattering_operator')then
			if(adjustl(trim(arg)).eq.'ptch') nscat= 0
			if(adjustl(trim(arg)).eq.'ener') nscat= 1
			if(adjustl(trim(arg)).eq.'both') nscat= 2
		endif
		if(trim(text).eq.'particle_type')then
			if(adjustl(trim(arg)).eq.'ion') cs=1
			if(adjustl(trim(arg)).eq.'ele') cs=-1
		endif
		if(trim(text).eq.'ions_mass') read(arg,*) mi
		if(trim(text).eq.'ions_charge') read(arg,*) zi
		if(trim(text).eq.'ions_temperature') read(arg,*) ntmpi
		if(trim(text).eq.'electrons_temperature') read(arg,*) ntmpe
		if(trim(text).eq.'temperature_profile')then
			if(adjustl(trim(arg)).eq.'flat') ntemp=0
			if(adjustl(trim(arg)).eq.'gaus') ntemp=1
		endif
		if(trim(text).eq.'density') read(arg,*) deno
		if(trim(text).eq.'density_profile')then
			if(adjustl(trim(arg)).eq.'flat') ndens=0
			if(adjustl(trim(arg)).eq.'para') ndens=1
		endif
		if(trim(text).eq.'mass') read(arg,*) mt
		if(trim(text).eq.'charge') read(arg,*) zt
		if(trim(text).eq.'pitch_angle') read(arg,*) ptchgrad
		if(trim(text).eq.'energy') read(arg,*) nengo
		if(trim(text).eq.'energy_distribution') then 
			if(adjustl(trim(arg)).eq.'mon') nends= 0
			if(adjustl(trim(arg)).eq.'max') nends= 1
		endif
		if(trim(text).eq.'radius') read(arg,*) nrado
		if(trim(text).eq.'theta') read(arg,*) thtgrad
		if(trim(text).eq.'zeta') read(arg,*) zetgrad
		if(trim(text).eq.'traj.plt') lngargtraj= line
		if(trim(text).eq.'poin.plt') lngargpoin= line
		if(trim(text).eq.'dist.plt') lngargdist= line
		if(trim(text).eq.'disp.plt') lngargdisp= line
		if(trim(text).eq.'coef.plt') lngargcoef= line
		if(eof.lt.0) exit
	enddo
	call userwrite
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine userwrite
!------------------------------------------!
! Read user output commands from data file !
!------------------------------------------!
!=========================================================!
	use shared,only:&
		pn01,&       !tim traj
		pn02,&       !rad traj
		pn03,&       !tht traj
		pn04,&       !zet traj
		pn05,&       !par traj
		pn06,&       !per traj
		pn07,&       !eng traj
		pn08,&       !pch traj
		pn09,&       !rcy traj
		pn10,&       !zcy traj
		pn11,&       !rad poin
		pn12,&       !tht poin
		pn13,&       !rcy poin
		pn14,&       !zcy poin
		pn21,&       !tim disp
		pn22,&       !dsr disp
		pn23,&       !dse disp
		pn24,&       !dsz disp
		pn25,&       !skr disp
		pn26,&       !ske disp
		pn27,&       !skz disp
		pn28,&       !flr disp
		pn29,&       !fle disp
		pn30,&       !flz disp
		pn31,&       !den coef
		pn32,&       !pot coef
		pn33,&       !npr coef
		pn34,&       !dfr coef
		pn35,&       !dfe coef
		pn36,&       !dfp coef
		pn37,&       !sdr coef
		pn38,&       !sde coef
		pn39,&       !sdp coef
		pn41,&       !rad dist
		pn42,&       !tht dist
		pn43,&       !zet dist
		trgtdum,&    !dummy target
		trgttraj,&   !traj targets
		trgtpoin,&   !poin targets
		trgtdist,&   !dist targets
		trgtdisp,&   !disp targets
		trgtcoef,&   !coef targets
		lngargtraj,& !traj.plt user string
		lngargpoin,& !poin.plt user string
		lngargdist,& !dist.plt user string
		lngargdisp,& !disp.plt user string
		lngargcoef   !coef.plt user string
	implicit none
	integer(8):: n   !index
	integer(8):: i   !index
	integer(8):: ind !index
	character(64):: str !string
!=========================================================!

	!traj.plt
	call strip(lngargtraj,str,n)
	allocate(trgttraj(n))
	pn01=> trgtdum
	pn02=> trgtdum
	pn03=> trgtdum
	pn04=> trgtdum
	pn05=> trgtdum
	pn06=> trgtdum
	pn07=> trgtdum
	pn08=> trgtdum
	pn09=> trgtdum
	pn10=> trgtdum
	ind= index(str,'tim')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn01=> trgttraj(i)
	ind= index(str,'rad')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn02=> trgttraj(i)
	ind= index(str,'tht')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn03=> trgttraj(i)
	ind= index(str,'zet')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn04=> trgttraj(i)
	ind= index(str,'par')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn05=> trgttraj(i)
	ind= index(str,'per')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn06=> trgttraj(i)
	ind= index(str,'eng')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn07=> trgttraj(i)
	ind= index(str,'pch')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn08=> trgttraj(i)
	ind= index(str,'rcy')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn09=> trgttraj(i)
	ind= index(str,'zcy')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn10=> trgttraj(i)

	!poin.plt
	call strip(lngargpoin,str,n)
	allocate(trgtpoin(n))
	pn11=> trgtdum
	pn12=> trgtdum
	pn13=> trgtdum
	pn14=> trgtdum
	ind= index(str,'rad')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn11=> trgtpoin(i)
	ind= index(str,'tht')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn12=> trgtpoin(i)
	ind= index(str,'rcy')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn13=> trgtpoin(i)
	ind= index(str,'zcy')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn14=> trgtpoin(i)

	!dist.plt
	call strip(lngargdist,str,n)
	allocate(trgtdist(n))
	pn41=> trgtdum
	pn42=> trgtdum
	pn43=> trgtdum
	ind= index(str,'rad')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn41=> trgtdist(i)
	ind= index(str,'eng')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn42=> trgtdist(i)
	ind= index(str,'pch')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn43=> trgtdist(i)

	!disp.plt
	call strip(lngargdisp,str,n)
	allocate(trgtdisp(n))
	pn21=> trgtdum
	pn22=> trgtdum
	pn23=> trgtdum
	pn24=> trgtdum
	pn25=> trgtdum
	pn26=> trgtdum
	pn27=> trgtdum
	pn28=> trgtdum
	pn29=> trgtdum
	pn30=> trgtdum
	ind= index(str,'tim')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn21=> trgtdisp(i)
	ind= index(str,'dsr')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn22=> trgtdisp(i)
	ind= index(str,'dse')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn23=> trgtdisp(i)
	ind= index(str,'dsz')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn24=> trgtdisp(i)
	ind= index(str,'skr')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn25=> trgtdisp(i)
	ind= index(str,'ske')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn26=> trgtdisp(i)
	ind= index(str,'skz')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn27=> trgtdisp(i)
	ind= index(str,'flr')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn28=> trgtdisp(i)	
	ind= index(str,'fle')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn29=> trgtdisp(i)
	ind= index(str,'flz')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn30=> trgtdisp(i)	
		
	!coef.plt
	call strip(lngargcoef,str,n)
	allocate(trgtcoef(n))
	pn31=> trgtdum
	pn32=> trgtdum
	pn33=> trgtdum
	pn34=> trgtdum
	pn35=> trgtdum
	pn36=> trgtdum
	ind= index(str,'den')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn31=> trgtcoef(i)
	ind= index(str,'pot')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn32=> trgtcoef(i)
	ind= index(str,'flx')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn33=> trgtcoef(i)
	ind= index(str,'dfr')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn34=> trgtcoef(i)
	ind= index(str,'dfe')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn35=> trgtcoef(i)
	ind= index(str,'dfp')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn36=> trgtcoef(i)
	ind= index(str,'sdr')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn37=> trgtcoef(i)
	ind= index(str,'sde')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn38=> trgtcoef(i)
	ind= index(str,'sdp')
	i= 1+ (ind -1)/3
	if(ind.gt.0) pn39=> trgtcoef(i)
	return
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine open_files(outdir)
!-----------------!
! Open data files !
!-----------------!
!=========================================================!
	use shared,only: &
		ntype,&      !input simulation swtich traj/poin/ense
		ndist,&      !input number of dist points
		lngargtraj,& !traj.plt user string
		lngargpoin,& !poin.plt user string
		lngargdisp,& !disp.plt user string
		lngargdist,& !dist.plt user string
		lngargcoef   !coef.plt user string
	implicit none
	logical:: ex
	integer(8):: j         !index
	integer(8):: n         !index
	character(64):: label  
	character(64):: narg
	character(64):: outdir !output directory
!=========================================================!
	
	write(label,*) trim(outdir)//'/data.log'
	open(36,file=adjustl(label),status='replace')
	if(ntype.eq.1)then	
		write(label,*) trim(outdir)//'/traj.plt'
		open(43,file=adjustl(label),status='replace')
		write(43,*) '# ',trim(lngargtraj)
	endif
	if(ntype.eq.2)then	
		write(label,*) trim(outdir)//'/disp.plt'
		open(50,file=adjustl(label),status='replace')
		write(50,*) '# ',trim(lngargdisp)
		write(label,*) trim(outdir)//'/rz.plt'
		do j=1,ndist
			write(narg,*) j
			label= trim(outdir)//'/dist_'//trim(adjustl(narg))//'.plt'
			n= 60 +j
			open(n,file=trim(label),status='replace')
		enddo
		write(label,*) trim(outdir)//'/coef.plt'
		inquire(file=adjustl(label),exist=ex)
		if(ex) open(70,file=adjustl(label),status='old',position='append')
		if(.not.ex)then
			open(70,file=adjustl(label),status='new',action='write')
			write(70,*) '# ',trim(lngargcoef)
		endif
	endif
	if(ntype.eq.3)then
		write(label,*) trim(outdir)//'/poin.plt'
		open(48,file=adjustl(trim(label)),status='replace')
		write(48,*) '# ',trim(lngargpoin)
	endif
	return	
end subroutine
!---------------------------------------------------------!


!---------------------------------------------------------!
subroutine info
!---------------!
! Displays info !
!---------------!
	implicit none
	write(0,'(a54)') ' '
	write(0,'(a54)') adjustl('GCAF= Guiding Center Analytic Field                   ')
	write(0,'(a54)') adjustl('This is a guiding center code developed by R. Saavedra')
	write(0,'(a54)') adjustl('Email: saavestro@gmail.com                            ')
	call rcrd_log(0)
end subroutine
!---------------------------------------------------------!
