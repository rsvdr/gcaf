!---------------------------------------------------------------------!
! GCAF= Guiding Center Analytic Field code by R. Saavedra             !
! Email of the author: saavestro@gmail.com                            !
! For more information please visit https://github.com/Saavestro/gcaf !
!---------------------------------------------------------------------!
!---------------------------------------------------------!
module shared
!--------------------------------!
! Shared constants and variables !
!--------------------------------!
	implicit none
	
	!CONSTANTS
	integer(8),parameter::    idm= 12000        !array dimension
	real(8),parameter::		  pi= 3.14159265d0  !pi number
	real(8),parameter::		  keverg= 1.6022d-9 !keV to erg
	real(8),parameter::		  btkg= 1.d+4       !10 kG
	real(8),parameter::		  rq= 4.8032d-10    !elec charge
	real(8),parameter::		  rc= 2.9979d+10    !speed of light
	real(8),parameter::		  mp= 1.6726d-24    !proton mass
	character(32),parameter:: fm='(18es12.4e2)' !output format

	!INPUT DATA VARIABLES
	integer(8):: 	ntype 	   !input simulation swtich traj/poin/ense
	integer(8)::	nconf  	   !input configuration switch tkmk/stell
	integer(8):: 	ncoll  	   !input collisions switch on/off
	integer(8):: 	npert 	   !input perturbation switch on/off
	integer(8):: 	npart      !input number of particles
	integer(8):: 	ntran  	   !input number of toroidal transits
	integer(8)::    nsdel      !input step delta
	integer(8):: 	nfixs 	   !input step fix switch on/off
	real(8)::    	deler 	   !input step fix en delta
	integer(8)::	ntraj 	   !input number of traj points
	integer(8)::	ndisp 	   !input number of disp points
	integer(8):: 	ndist      !input number of dist points
	integer(8)::  	bm    	   !input on axis b field
	integer(8)::   	rmaj  	   !input major radius
	integer(8)::    rmin 	   !input minor radius
	real(8)::       namx  	   !input perturbation amplitude
	integer(8)::    nfrq 	   !input perturbation frequency
	integer(8)::    mmod 	   !input m mode number
	integer(8)::    nmod 	   !input n mode number
	integer(8)::    nphs 	   !input perturbation phase
	real(8)::    	npto       !input pt in keV
	integer(8)::    nepot      !input pt profile
	integer(8):: 	nspec 	   !scatter species ions/elec/both
	integer(8):: 	nscat  	   !scatter pitch/en/both
	integer(8)::    mi    	   !ions mass in mp
	integer(8)::	zi    	   !ions charge in protons charge
	integer(8)::   	ntmpi 	   !ions temp in keV
	integer(8)::	ntmpe 	   !elec temp in keV
	integer(8)::    ntemp      !temp profile flat/gaus
	real(8)::    	deno   	   !plasma dens in cm^-3
	integer(8)::    ndens      !dens profile flat/para
	integer(8)::	mt    	   !test mass in mp
	integer(8)::	zt    	   !test charge in protons charge
	integer(8)::    ptchgrad   !input pitch in grades
	real(8)::       nengo      !input initial en in keV		
	integer(8):: 	nends 	   !input en distribution mon/max
	integer(8)::    nrado  	   !input initial rad in cm
	integer(8)::    thtgrad    !theta in grades
	integer(8)::    zetgrad    !zeta in grades
	integer(8)::    nscrn  	   !screen output on/off
	integer(8)::    cs         !test type ions/elec	
	character(64)::	lngargtraj !traj.plt user string
	character(64)::	lngargpoin !poin.plt user string
	character(64)::	lngargdisp !disp.plt user string
	character(64)::	lngargdist !dist.plt user string
	character(64)::	lngargcoef !coef.plt user string			

	!SIMULATION
	real(8)::    					   trun  !total simulation time
	real(8)::    					   omego !characteristic time freq
	real(8)::						   leno  !characteristic length
	real(8)::   					   enon  !normalized initial en
	real(8)::    					   enaxo !on axis en 
	real(8)::						   dto   !initial time step
	real(8),dimension(:),allocatable:: dt    !time step
	real(8),dimension(:),allocatable:: time  !time variable

    !FIELD
   	real(8)::       				   prtfrq !perturbation frequency
	real(8)::						   pto    !pt magnitude
	real(8),dimension(:),allocatable:: bf     !b field 
	real(8),dimension(:),allocatable:: dbdp   !b pol derivative
	real(8),dimension(:),allocatable:: dbdt   !b tht derivative
	real(8),dimension(:),allocatable:: dbdz   !b zet derivative
	real(8),dimension(:),allocatable:: al     !alfa perturbation 
	real(8),dimension(:),allocatable:: dadp   !al pol derivative
	real(8),dimension(:),allocatable:: dadt   !al tht derivative
	real(8),dimension(:),allocatable:: dadz   !al zet derivative
	real(8),dimension(:),allocatable:: dadtim !al zet derivative
	real(8),dimension(:),allocatable:: pt     !electrostatic potential
	real(8),dimension(:),allocatable:: dptdp  !pt pol derivative
	real(8),dimension(:),allocatable:: dptdt  !pt tht derivative
	real(8),dimension(:),allocatable:: dptdz  !pt zet derivative
	real(8),dimension(:),allocatable:: cg     !poloidal current
	real(8),dimension(:),allocatable:: cgp    !cg pol derivative
	real(8),dimension(:),allocatable:: ci     !toroidal current
	real(8),dimension(:),allocatable:: cip    !ci pol derivative
	real(8),dimension(:),allocatable:: qf     !q function

    !TEST PARTICLE
	real(8)::		                   ptcho  !initial pitch in cos
	real(8)::						   eno    !initial en in keV
	real(8)::					       polo   !initial poloidal flux
	real(8)::						   theto  !initial theta
	real(8)::						   zeto   !initial zeta
	real(8)::                          polw   !on wall poloidal flux
	real(8),dimension(:),allocatable:: ptch   !pitch
	real(8),dimension(:),allocatable:: en     !energy
	real(8),dimension(:),allocatable:: pol    !poloidal flux
	real(8),dimension(:),allocatable:: tht    !theta
	real(8),dimension(:),allocatable:: zet    !zeta
	real(8),dimension(:),allocatable:: rho    !parallel gyroradius
	real(8),dimension(:),allocatable:: mu     !b moment
	real(8),dimension(:),allocatable:: pz     !zeta conjugate momenta
	real(8),dimension(:),allocatable:: toten  !total en
	real(8),dimension(:),allocatable:: totenn !total en normalized

    !ENSEMBLE
	real(8),dimension(:),allocatable:: dft       !time for disp
	real(8),dimension(:),allocatable:: rav       !rad av
	real(8),dimension(:),allocatable:: rsqav     !rad sq av
	real(8),dimension(:),allocatable:: rcuav     !rad cu av
	real(8),dimension(:),allocatable:: rfoav     !rad fo av
	real(8),dimension(:),allocatable:: enav      !en av
	real(8),dimension(:),allocatable:: ensqav    !en sq av
	real(8),dimension(:),allocatable:: encuav    !en cu av
	real(8),dimension(:),allocatable:: enfoav    !en fo av
	real(8),dimension(:),allocatable:: pzav      !pz av
	real(8),dimension(:),allocatable:: pzsqav    !pz sq av
	real(8),dimension(:),allocatable:: pzcuav    !pz cu av
	real(8),dimension(:),allocatable:: pzfoav    !pz fo av
	real(8),dimension(:),allocatable:: rone      !rad one part
	real(8),dimension(:),allocatable:: rsqone    !rad sq one part
	real(8),dimension(:),allocatable:: rcuone    !rad cu one part
	real(8),dimension(:),allocatable:: rfoone    !rad fo one part
	real(8),dimension(:),allocatable:: enone     !en one part
	real(8),dimension(:),allocatable:: ensqone   !en sq one part
	real(8),dimension(:),allocatable:: encuone   !en cu one part
	real(8),dimension(:),allocatable:: enfoone   !en fo one part
	real(8),dimension(:),allocatable:: pzone     !pz one part 
	real(8),dimension(:),allocatable:: pzsqone   !pz sq one part
	real(8),dimension(:),allocatable:: pzcuone   !pz cu one part
	real(8),dimension(:),allocatable:: pzfoone   !pz fo one part
	real(8),dimension(:),allocatable:: rprv      !rad previous

    !RECORDS
    integer(8)::                         flux     !particle flux
	real(8)::                            r1       !poin lower bound
	real(8)::                            r2       !poin upper bound
	real(8)::                            zplo     !poin zeta plane
	real(8),dimension(:),allocatable::   pold     !old pol for poin
	real(8),dimension(:),allocatable::   told     !old tht for poin
	real(8),dimension(:),allocatable::   zold     !old zet for poin
	real(8),dimension(:),allocatable::   distt    !time dist
	real(8),dimension(:,:),allocatable:: distr    !rad dist
	real(8),dimension(:,:),allocatable:: diste    !en dist
	real(8),dimension(:,:),allocatable:: distp    !ptch dist
	real(8),dimension(:,:),allocatable:: rcrdtraj !record traj points
	real(8),dimension(:,:),allocatable:: rcrdpoin !record poin points
	real(8),pointer:: pn01 !tim traj
	real(8),pointer:: pn02 !rad traj
	real(8),pointer:: pn03 !tht traj
	real(8),pointer:: pn04 !zet traj
	real(8),pointer:: pn05 !par traj
	real(8),pointer:: pn06 !per traj
	real(8),pointer:: pn07 !eng traj
	real(8),pointer:: pn08 !pch traj
	real(8),pointer:: pn09 !rcy traj
	real(8),pointer:: pn10 !zcy traj
	real(8),pointer:: pn11 !rad poin
	real(8),pointer:: pn12 !tht poin
	real(8),pointer:: pn13 !rcy poin
	real(8),pointer:: pn14 !zcy poin
	real(8),pointer:: pn15 
	real(8),pointer:: pn16
	real(8),pointer:: pn17
	real(8),pointer:: pn18
	real(8),pointer:: pn19
	real(8),pointer:: pn20 
	real(8),pointer:: pn21 !tim disp
	real(8),pointer:: pn22 !dsr disp
	real(8),pointer:: pn23 !dse disp
	real(8),pointer:: pn24 !dsz disp
	real(8),pointer:: pn25 !skr disp
	real(8),pointer:: pn26 !ske disp
	real(8),pointer:: pn27 !skz disp
	real(8),pointer:: pn28 !flr disp
	real(8),pointer:: pn29 !fle disp
	real(8),pointer:: pn30 !flz disp
	real(8),pointer:: pn31 !den coef
	real(8),pointer:: pn32 !pot coef
	real(8),pointer:: pn33 !flx coef
	real(8),pointer:: pn34 !dfr coef
	real(8),pointer:: pn35 !dfe coef
	real(8),pointer:: pn36 !dfp coef
	real(8),pointer:: pn37 !sdr coef
	real(8),pointer:: pn38 !sde coef
	real(8),pointer:: pn39 !sdp coef
	real(8),pointer:: pn40 
	real(8),pointer:: pn41 !rad dist
	real(8),pointer:: pn42 !tht dist
	real(8),pointer:: pn43 !zet dist
	real(8),pointer:: pn44
	real(8),target::                          trgtdum  !dummy target
	real(8),target,dimension(:),allocatable:: trgttraj !traj targets
	real(8),target,dimension(:),allocatable:: trgtpoin !poin targets
	real(8),target,dimension(:),allocatable:: trgtdist !dist targets
	real(8),target,dimension(:),allocatable:: trgtdisp !disp targets
	real(8),target,dimension(:),allocatable:: trgtcoef !coef targets
    character(4),dimension(:),allocatable:: prnt
end module
!---------------------------------------------------------!
