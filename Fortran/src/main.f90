!---------------------------------------------------------------------!
! GCAF= Guiding Center Analytic Field code by R. Saavedra             !
! Email of the author: saavestro@gmail.com                            !
! For more information please visit https://github.com/Saavestro/gcaf !
!---------------------------------------------------------------------!
program main
!=========================================================!
	use shared,only:&
		ntype,& !simulation type traj/poin/ense
		ncoll,& !collisions on/off
		npart,& !number of particles
		nscrn,& !screen output on/off
		nfixs,& !step fix on/off
		time,&  !time variable
		trun,&  !total simulation time
		ntraj   !number of traj record points
	implicit none
	integer(8):: k     		!particle index
	integer(8):: nskip 		!trajectory record skip
	integer(8):: mskip 		!dispersion record skip
	integer(8):: nend       !end simulation
	integer(8):: nthreads   !number of OMP threads
	integer(8):: start_time !start wall time
	integer(8):: end_time   !end wall time
	integer(8):: clock_rate !wall time variable
!=========================================================!
	!$ nthreads=4
	!$ call omp_set_num_threads(nthreads)

	!INITIALIZE SIMULATION
	call cmdlin
	call system_clock(count=start_time)
	call inivec
	call inivar
	do k=1,npart
		call set(k)
	enddo
	
	if(ntype.ne.1) go to 22
	!TRAJECTORY LOOP
	k= 1
	nend= 0
	nskip= 1
	mskip= 1
	call rcrd_traj(k,nskip)
	do while(nend.ne.1)
		call onestep(k)
		if(nfixs.eq.1) call stepfix(k)
		call update(k,nskip,mskip,nend)
		if(ncoll.eq.1) call scatr(k)
		if(time(k).gt.nskip*trun/ntraj) call rcrd_traj(k,nskip)
	enddo
	go to 23
	
	22 continue
	!ENSEMBLE LOOP
	!$omp parallel do private(k,nskip,mskip,nend)
	do k=1,npart
		11 continue
		nend= 0
		nskip= 1
		mskip= 1
		do while(nend.ne.1)
			call onestep(k)
			if(nfixs.eq.1) call stepfix(k)
			call update(k,nskip,mskip,nend)
			if(ncoll.eq.1) call scatr(k)
			if(ntype.eq.3) call write_poin(k)
		enddo
		call update_ensemble
		if(nscrn.eq.1) write(0,'(a9,i4,a9)') 'particle ',k,' finished'
	enddo
	!$omp end parallel do
	go to 23
	
	23 continue
	!RECORDS
	if(ntype.eq.1) call write_traj
	if(ntype.eq.2) then
		call write_disp	
		call write_dist
	endif

	!END SIMULATION
	call finvec
	call system_clock(count=end_time,count_rate=clock_rate)
	if (clock_rate.gt.0) write(0,'(A10,F10.2)') 'CPU TIME= ', &
		dble(end_time -start_time)/dble(clock_rate)
end program
