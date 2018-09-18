!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018-08-25
!--- Acknowledgement to Liu group, of PKU
!*

!-- ******************************************************************************************
!--                      Mapping Schrodinger Equation to quaternion space
!-- ******************************************************************************************
!-- Basic definitions:
!------ quaternion units: { 1, i, j, k }
!------ >> 1^2 = 1, i^2 = j^2 = k^2 = -1
!------ >> ij=-ji=k, jk=-kj=i, ki=-ik=j
!-- Operations:
!------ conjugattion : x=a+bi+cj+dk  ==>  x*=a-bi-cj-dk
!------ normsquare   : x=a+bi+cj+dk  ==>  |x|^2 = x*x = a^2+b^2+c^2+d^2
!-- Analogy to SE:
!------ the Schrodinger Equation is the section (complex number space) of the quaternion space.
!------ let:
!--------- u = cos(theta)*i + sin(theta)cos(phi) * j + sin(theta)sin(phi) * k
!--------- and n = [ cos(theta), sin(theta)cos(phi), sin(theta)sin(phi) ], to be direction cosine
!------ then:
!--------- u^2 = -1, so u is the general imaginary unit in quaternion space
!------ and Schrodinger euqation re-formulate as:
!--------- u*dc/dt = H c
!------ where u is the general "imaginary" unit in quaternion space, and c is of quaternion.
!-- ******************************************************************************************

!-- Some implementary
!------ 1. Models: ModelI~ModelVI, same as Liu's paper. Model VII, is linear mixture of propagator,
!------    while ModelVIII, trys to combine different mapping methods (pseudo-rotation). and 
!------    Model IX is Spin-mapping model, finally ModelX is the Quarternion version model.

module map4
implicit none
	!-- 1) HAMILTONIAN SETTINGS
	integer :: map_mtxsize
	!-- allocatable objects of H-matrix etc.
	real(8), dimension(:,:), allocatable :: map_mtxH ! Hamiltonian
	real(8), dimension(:,:), allocatable :: map_mtxV ! kinds of valuables
	real(8), dimension(:,:), allocatable :: map_mtxO ! matixe of e^(-H)
	real(8), dimension(:), allocatable   :: map_mtxN   ! Populations
	!-- modification arguments of Hamiltonian, if not set the value < 0
	real(8) :: map_Hnondiag = -1.0
	
	!-- 2) MODEL ARGUMENTS
	!-- outer arguments for different models
	integer :: map_model
	!-- inner arguments of different models
	integer :: map_vardim ! dimension of varibles
	integer :: map_initmd ! initialization method
	integer :: map_analmd ! analyzation method
	
	!-- 3) DYNAMICS SETTINGS
	!-- if imaginary, set map_imag = 1, otherwise set 0 for real dynamics
	integer :: map_imag = 0
	!-- general for dynamcs simulation
	real(8) :: map_dtime
	integer :: map_nstep
	integer :: map_npass
	integer :: map_nsamp
	
	!-- 4) INITIALIZATION SNETTINGS
	!-- the initial populattion distribution, if map_istat == 0, then arrange the pupolation by default
	!----- otherwise, populate on i-th state
	integer :: map_istat
	!-- initial phase angle
	real(8) :: map_angle
	
	!-- 5) OTHERS
	!-- quaternion model unit -- direction cosine
	real(8), dimension(3) :: map_dc
	!-- thermodynamics arguments
	real(8) :: map_temp
	real(8) :: map_beta
	!-- [math]
	real(8) :: map_pi  = 3.14159265358979323846
	real(8) :: map_2pi = 2.0_8 * 3.14159265358979323846
	!-- temporary varibles
	real(8) :: map_cos
	real(8) :: map_sin
	

contains
	!-- defination of various dot operation
	!------ note that: matrix H and vector x should have same size in one direction, here we don't verify anymore
	!------ note that: H_dot = H_ndot + H_ddot     :     ( nondiagonal + diagonal )
	real(8) function H_dot(n,flag)
		integer, intent(in) :: n, flag
		H_dot = dot_product( map_mtxH(n,:), map_mtxV(flag,:) )
	end function H_dot
	
	real(8) function H_ndot(n,flag)
		integer, intent(in) :: n, flag
		H_ndot = dot_product( map_mtxH(n,:), map_mtxV(flag,:) )
		H_ndot = H_ndot - map_mtxH(n,n)*map_mtxV(flag,n)
	end function H_ndot
	
	real(8) function H_ddot(n,flag)
		integer, intent(in) :: n, flag
		H_ddot = map_mtxH(n,n)*map_mtxV(flag,n)
	end function H_ddot
	
	real(8) function H_xdot(n,flag) ! half-diagonal dot operation
		integer, intent(in) :: n, flag
		H_xdot = dot_product( map_mtxH(n,:), map_mtxV(flag,:) )
		H_xdot = H_xdot - 0.5_8 * map_mtxH(n,n)*map_mtxV(flag,n)
	end function H_xdot
	
	!-- spin compoment in x-y plane, of n-th sub-space
	real(8) function spinm(n)
	    integer, intent(in) :: n
	    real(8) :: S2
	    S2 = 0.75
	    spinm = dsqrt( S2 - map_mtxV(2,n)**2 )
	end function spinm
	
	!-- cos(thetai - thetaj)
	real(8) function cosd(i,j)
	    integer, intent(in) :: i,j
	    cosd = dcos( map_mtxV(1,i) - map_mtxV(1,j) )
	end function cosd
	
	!-- sin(thetai - thetaj)
	real(8) function sind(i,j)
	    integer, intent(in) :: i,j
	    sind = dsin( map_mtxV(1,i) - map_mtxV(1,j) )
	end function sind
	
	!-- calculate population from i-state to j-state
	real(8) function map_pop(popi,popf)
	    integer, intent(in) :: popi,popf
	    integer :: i
	    map_pop = 0
	    do i=popi, popf
	        if(map_analmd.eq. 1) then
	            map_pop = map_pop + map_mtxV(1,i)*map_mtxV(4,i) - map_mtxV(3,i)*map_mtxV(2,i)
    	    elseif(map_analmd .eq. 2) then
    	        map_pop = map_pop + 0.5_8 * (map_mtxV(1,i)**2 + map_mtxV(2,i)**2)
    	    elseif(map_analmd .eq. 3) then
    	        map_pop = map_pop + 0.25_8 * ( (map_mtxV(1,i)+map_mtxV(4,i))**2 + (map_mtxV(3,i)-map_mtxV(2,i))**2 )
    	    else
    	        print *,"map_pop() not support"
    	        stop
    	    end if
	    enddo
	end function map_pop
	
	!-- renormalization
	subroutine map_renorm()
		real(8) :: pop
		if(map_imag .eq. 1) then
			if(map_analmd .eq. 2) then
				pop = 0
			    do i=1,map_mtxsize
			        pop = pop + 0.5*(map_mtxV(1,i)**2 + map_mtxV(2,i)**2)
			    end do
			    map_mtxV(:,:) = map_mtxV(:,:)/sqrt(pop)
			end if
		end if
	end subroutine map_renorm
	
	!-- if imaginary dynamics, to obtain excited components, should eliminate the ground state components
	subroutine map_eraser(elimfile)
	    character(*), intent(in) :: elimfile
	    real(8), dimension(map_vardim,map_mtxsize) :: tempvs
	    real(8) :: prj_re, prj_im
	    integer :: i, myiostat
	    if(map_imag.eq. 1) then
	        !-- read file for elimination
		    open(unit=23,file=trim(elimfile),status="old")
		    do i=1,map_mtxsize
		        read(23,*,iostat=myiostat) tempvs(:,i)
		        if(myiostat<0) exit
		    enddo
		    close(unit=23)
		    
		    !-- only analmd2 supported
            if(map_analmd .eq. 2) then
                !-- calculate porjection of VS_excited and VS_ground
			    prj_re = 0
			    prj_im = 0
			    do i=1,map_mtxsize
			        prj_re = prj_re + tempvs(1,i)*map_mtxV(1,i) + tempvs(2,i)*map_mtxV(2,i)
			        prj_im = prj_im - tempvs(2,i)*map_mtxV(1,i) + tempvs(1,i)*map_mtxV(2,i)
			    enddo
			    
			    !-- do elimination
			    do i=1,map_mtxsize
			        map_mtxV(1,i) = map_mtxV(1,i) - (prj_re/2.0) * tempvs(1,i)
			        map_mtxV(2,i) = map_mtxV(2,i) - (prj_im/2.0) * tempvs(2,i)
			    enddo
			end if
	    end if
	    return
	end subroutine map_eraser
	
	
	!-- random seed
	subroutine init_seed()
	    integer :: n, ival(8), v(3), i
	    integer, allocatable :: seed(:)
	    call date_and_time(values=ival)
	    v(1) = ival(8) + 2048*ival(7)
	    v(2) = ival(6) + 64*ival(5)     ! value(4) isn't real(kind=8)ly 'random'
	    v(3) = ival(3) + 32*ival(2) + 32*8*ival(1)
     	call random_seed(size=n)
	    allocate(seed(n))
	    call random_seed()   ! Give the seed an implementation-dependent kick
	    call random_seed(get=seed)
	    do i=1, n
        	seed(i) = seed(i) + v(mod(i-1, 3) + 1)
      	enddo
      	call random_seed(put=seed)
      	deallocate(seed)
    end subroutine init_seed


	!-- read from 'map.rc' and 'Hsave.dat' to initialize arguments list and Hamiltonian matirx
	subroutine map_initmap()
	implicit none
		character(len=20), dimension(2) :: pairs
		integer :: i,j
		integer :: my_count, my_iostat
		logical :: my_exist
		real(8) :: tmp_norm
		
		!-- read arguments from 'map.rc' file
		!------ note that: my_count indicates the number of valid initialized object
		inquire( file=trim('map.rc'), exist=my_exist )
		if (my_exist .eqv. .true.) then
		    open( unit=10, file=trim('map.rc') )
		    my_count = 0
		    do while (.true.)
		        read(10,*,iostat=my_iostat) pairs
		        if (my_iostat < 0) exit
		        my_count = my_count + 1
		        select case (pairs(1))
		            case ('mapmodel')
		                read(pairs(2),*) map_model
		                select case (map_model)
		                    case (1)
		                        map_vardim = 4
		                        map_initmd = 1
		                        map_analmd = 1
		                    case (2)
		                        map_vardim = 2
		                        map_initmd = 2
		                        map_analmd = 2
		                    case (3)
		                        map_vardim = 4
		                        map_initmd = 1
		                        map_analmd = 3
		                    case (4)
		                        map_vardim = 4
		                        map_initmd = 1
		                        map_analmd = 3
		                    case (5)
		                        map_vardim = 4
		                        map_initmd = 1
		                        map_analmd = 3
		                    case (6)
		                        map_vardim = 4
		                        map_initmd = 1
		                        map_analmd = 1
		                    case (7)
		                        map_vardim = 4
		                        map_initmd = 1
		                        map_analmd = 3
		                    case (8)
		                        map_vardim = 4
		                        map_initmd = 3
		                        map_analmd = 4
		                    case (9)
		                        map_vardim = 2
		                        map_initmd = 4
		                        map_analmd = 5
		                    case (10)
		                        map_vardim = 4
		                        map_initmd = 5
		                        map_analmd = 6
		                    case default
		                        print *, "init failed of model"
		                        stop
		                end select
		            case ('phsangle')
		            	read(pairs(2),*) map_angle
		            case ('mapdtime')
		            	read(pairs(2),*) map_dtime
		            case ('mapbeta')
		            	read(pairs(2),*) map_beta
		            	map_temp = 1.0/map_beta
		            case ('mapnstep')
		            	read(pairs(2),*) map_nstep
		            case ('mapnsamp')
		            	read(pairs(2),*) map_nsamp
		            case ('inistate')
		            	read(pairs(2),*) map_istat
		            case ('hnondiag')
		            	read(pairs(2),*) map_Hnondiag
		            case ('mapimag')
		            	read(pairs(2),*) map_imag
		            case ('dc1')
		                read(pairs(2),*) map_dc(1)
		            case ('dc2')
		                read(pairs(2),*) map_dc(2)
		            case ('dc3')
		                read(pairs(2),*) map_dc(3)    
		            case default
		            	print *, "warning, arguments in 'map.rc' doesn't match!"
		            	my_count = my_count - 1
		        end select
		    end do
		else
			print *, "error, 'map.rc' doesn't exist!"
			stop
		end if
		
		tmp_norm = sqrt( map_dc(1)**2 + map_dc(2)**2 + map_dc(3)**2 )
		map_dc(1) = map_dc(1)/tmp_norm
		map_dc(2) = map_dc(2)/tmp_norm
		map_dc(3) = map_dc(3)/tmp_norm
		print *, "direction cosine: ", map_dc
		
		!-- verify the completeness of arguments
		if(my_count .ge. 12) then
			print *, "initialized from 'map.rc' successfully!"
		else
			print *, "error, initialization from 'map.rc' fails!"
		    stop
		end if
		
		!-- read the information of Hamitonian of a given state-space
		!------ 1) inquire existing of the Hamiltonian
		inquire(file=trim('Hsave.dat'), exist=my_exist)
        if (my_exist .eqv. .false.) then
           	print *, "error, Hsave.dat doesn't exist"
           	stop
        end if
        open(unit=10, file=trim('Hsave.dat'), status='old')
        
        !------ 2) initialize the size of state space
        read(10, *) map_mtxsize
        allocate(map_mtxH(map_mtxsize, map_mtxsize))
        allocate(map_mtxO(map_mtxsize, map_mtxsize))
        allocate(map_mtxV( map_vardim , map_mtxsize))   ! alternatively, just allocate size of (4, map_mtxzize)
        allocate(map_mtxN(map_mtxsize))
        
        !------ 3) read Hamitonian datas
        !------------- note that: my_count indicates the line-number of hamiltonian
        do my_count=1,map_mtxsize
        	read(10,*,iostat=my_iostat) map_mtxH(my_count,:)
		    if (my_iostat < 0) exit
        end do
        !-- modification of Hamiltonian
        if(map_Hnondiag .ge. 0) then
        	do i=1,map_mtxsize
        		do j=i+1,map_mtxsize
        			map_mtxH(i,j) = map_Hnondiag
        			map_mtxH(j,i) = map_Hnondiag
        		end do
        	end do
        end if
        close(unit=10)
        
        !-- show and verify the Hamiltonian
        do my_count=1,map_mtxsize
        	print *, map_mtxH(my_count,:)
        end do
        
        !-- initial the O-matrix at infinite temperature
        !-- ??
        map_mtxO(:,:) = 0.0
        do i=1,map_mtxsize
        	map_mtxO(i,i) = 1.0/map_mtxsize
        enddo
	end subroutine map_initmap
	
	
	
	!-- propagate the assistant valuables
	subroutine map_propagator(dtime)
	implicit none
		real(8), intent(in) :: dtime
		integer :: i,j
		real(8) :: irand, x, y, c, s, cc, ss, sc
		real(8) :: Elevel1, Elevel2
		
		!-- differet model, different propagator
		!------ note that:
		!------ 1) map_mtxV(1,:)   stands for x
		!------ 1) map_mtxV(2,:)   stands for Px
		!------ 1) map_mtxV(3,:)   stands for y
		!------ 1) map_mtxV(4,:)   stands for Py
		
		!-- REAL dynamics
		if(map_imag .eq. 0) then
		    select case (map_model)
			    case (1)
				    do i=1, map_mtxsize
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/2) * H_dot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/2) * H_dot(i,4)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime)   * H_dot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime)   * H_dot(i,2)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/2) * H_dot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/2) * H_dot(i,4)
				    end do
			    case (2)
			        do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/2) * H_dot(i,1)
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime)   * H_dot(i,2)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/2) * H_dot(i,1)
				    end do
			    case (3)
				    do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_ddot(i,1)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * H_xdot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_xdot(i,4)
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * H_ddot(i,2)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * H_xdot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_xdot(i,4)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_ddot(i,1)
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_ddot(i,3)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * H_xdot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/2) * H_xdot(i,2)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * H_ddot(i,4)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * H_xdot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/2) * H_xdot(i,2)
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_ddot(i,3)
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_ddot(i,1)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * H_xdot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_xdot(i,4)
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * H_ddot(i,2)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * H_xdot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_xdot(i,4)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_ddot(i,1)
				    end do
			    case (4)
				    do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_ddot(i,4)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_xdot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_xdot(i,3)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * H_ddot(i,2)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_xdot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_xdot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_ddot(i,4)
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * H_ddot(i,3)
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/2) * H_xdot(i,2)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * H_xdot(i,4)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * H_ddot(i,1)
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/2) * H_xdot(i,2)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * H_xdot(i,4)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * H_ddot(i,3)
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_ddot(i,4)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_xdot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_xdot(i,3)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * H_ddot(i,2)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_xdot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_xdot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_ddot(i,4)
				    end do
			    case (5)
				    do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_dot(i,4)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_dot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * H_dot(i,3)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * H_dot(i,2)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_dot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * H_dot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_dot(i,4)
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * H_dot(i,3)
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * H_dot(i,2)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/4) * H_dot(i,4)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * H_dot(i,1)
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * H_dot(i,2)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/4) * H_dot(i,4)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * H_dot(i,3)
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_dot(i,4)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_dot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * H_dot(i,3)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * H_dot(i,2)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_dot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * H_dot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * H_dot(i,4)
				    end do
			    case (6)
				    do i=1, map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_ddot(i,4)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_ndot(i,1) 
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_ndot(i,3)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/2) * H_ddot(i,2) 
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_ndot(i,1) 
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_ndot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_ddot(i,4)
					
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/2) * H_ddot(i,3) 
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/2) * H_ndot(i,2) 
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * H_ndot(i,4)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/1) * H_ddot(i,1)
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/2) * H_ndot(i,2)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * H_ndot(i,4)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/2) * H_ddot(i,3)
					
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_ddot(i,4)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_ndot(i,1) 
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_ndot(i,3)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/2) * H_ddot(i,2) 
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_ndot(i,1) 
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/4) * H_ndot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/4) * H_ddot(i,4)
				    end do
			    case (7)
			        irand = 0.5
				    do i=1,map_mtxsize
				        !call random_number(irand)
				        x = 2*irand
				        y = 2 - x
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*H_ndot(i,3) - y*H_ndot(i,2) )
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * (-x*H_ddot(i,3) + y*H_ddot(i,2) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*H_ndot(i,3) - y*H_ndot(i,2) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * ( x*H_ddot(i,3) - y*H_ddot(i,2) )
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * ( y*H_ndot(i,2) - x*H_ndot(i,3) )
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/4) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * ( x*H_ndot(i,2) - y*H_ndot(i,3) )
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/4) * ( y*H_ndot(i,4) + x*H_ndot(i,1) )
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * ( x*H_ddot(i,3) - y*H_ddot(i,2) )
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*H_ndot(i,3) - y*H_ndot(i,2) )
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * (-x*H_ddot(i,3) + y*H_ddot(i,2) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*H_ndot(i,3) - y*H_ndot(i,2) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
				    end do
				!-- WITH BUGS of case 8 --- model of mixed mapping
		        case (8)
		            c = sqrt(2.0)*map_cos
		            s = sqrt(2.0)*map_sin
		            cc = sqrt(2.0)*map_cos*map_sin
		            ss = sqrt(2.0)*map_sin*map_sin
		            sc = 1.0
				    call random_number(irand)
				    x = 1.0!x = 2*irand
				    y = 1.0!y = 2 - x
				    do i=1,map_mtxsize
				        map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*cc*H_ddot(i,1) + y*sc*H_ddot(i,4) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*cc*H_ndot(i,1) + y*sc*H_ndot(i,4) )
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*cc*H_ndot(i,3) - y*sc*H_ndot(i,2) )
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * (-x*cc*H_ddot(i,3) + y*sc*H_ddot(i,2) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*cc*H_ndot(i,1) + y*sc*H_ndot(i,4) )
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*cc*H_ndot(i,3) - y*sc*H_ndot(i,2) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*cc*H_ddot(i,1) + y*sc*H_ddot(i,4) )
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * ( x*sc*H_ddot(i,3) - y*ss*H_ddot(i,2) )
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * ( y*ss*H_ndot(i,2) - x*sc*H_ndot(i,3) )
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/4) * ( x*sc*H_ndot(i,1) + y*ss*H_ndot(i,4) )
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * ( x*sc*H_ddot(i,1) + y*ss*H_ddot(i,4) )
					    map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * ( x*sc*H_ndot(i,2) - y*ss*H_ndot(i,3) )
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime/4) * ( y*ss*H_ndot(i,4) + x*sc*H_ndot(i,1) )
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * ( x*sc*H_ddot(i,3) - y*ss*H_ddot(i,2) )
				    !end do
				    !do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*cc*H_ddot(i,1) + y*sc*H_ddot(i,4) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*cc*H_ndot(i,1) + y*sc*H_ndot(i,4) )
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*cc*H_ndot(i,3) - y*sc*H_ndot(i,2) )
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * (-x*cc*H_ddot(i,3) + y*sc*H_ddot(i,2) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*cc*H_ndot(i,1) + y*sc*H_ndot(i,4) )
					    map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*cc*H_ndot(i,3) - y*sc*H_ndot(i,2) )
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*cc*H_ddot(i,1) + y*sc*H_ddot(i,4) )
				    end do
			    !-- for spin mapping model
			    !-- map_mtxV(1,i):    q(i)  angular action
			    !-- map_mtxV(2,i):    m(i)  angular momenta
			    case (9)
				    do i=1, map_mtxsize
				        map_mtxV(1,i) = map_mtxV(1,i) + (dtime) * (1.0 + map_mtxV(2,i)) * map_mtxH(i,i)
				        do j=1, map_mtxsize
				            if(i.eq.j) cycle
				            map_mtxV(1,i) = map_mtxV(1,i) - (dtime) * (2.0*map_mtxV(2,i)/spinm(i)) * spinm(j) * cosd(i,j) *map_mtxH(i,j)
				        end do
				    end do
				    do i=1, map_mtxsize
				        do j=1, map_mtxsize
				            if(i.eq.j) cycle
				            map_mtxV(2,i) = map_mtxV(2,i) + (dtime) * spinm(i) * spinm(j) * sind(i,j) *map_mtxH(i,j)
				        end do
				    end do
				!-- quarternion mapping
				!------ performance: accuracy at dt=0.001 for three-state model
				case (10)
				    do i=1, map_mtxsize
				        map_mtxV(1,i) = map_mtxV(1,i) + (map_dtime/2.0) * &
				        	( map_dc(1)*H_dot(i,2) + map_dc(2)*H_dot(i,3) + map_dc(3)*H_dot(i,4) )
				        map_mtxV(2,i) = map_mtxV(2,i) + (map_dtime/4.0) * &
				        	(-map_dc(1)*H_dot(i,1) + map_dc(3)*H_dot(i,3) - map_dc(2)*H_dot(i,4) )
				        map_mtxV(4,i) = map_mtxV(4,i) + (map_dtime/2.0) * &
				        	(-map_dc(3)*H_dot(i,1) + map_dc(2)*H_dot(i,2) - map_dc(1)*H_dot(i,3) )
				        map_mtxV(2,i) = map_mtxV(2,i) + (map_dtime/4.0) * &
				        	(-map_dc(1)*H_dot(i,1) + map_dc(3)*H_dot(i,3) - map_dc(2)*H_dot(i,4) )
				        map_mtxV(3,i) = map_mtxV(3,i) + (map_dtime) * &
				        	(-map_dc(2)*H_dot(i,1) - map_dc(3)*H_dot(i,2) + map_dc(1)*H_dot(i,4) )
				        map_mtxV(2,i) = map_mtxV(2,i) + (map_dtime/4.0) * &
				        	(-map_dc(1)*H_dot(i,1) + map_dc(3)*H_dot(i,3) - map_dc(2)*H_dot(i,4) )
				        map_mtxV(4,i) = map_mtxV(4,i) + (map_dtime/2.0) * &
				        	(-map_dc(3)*H_dot(i,1) + map_dc(2)*H_dot(i,2) - map_dc(1)*H_dot(i,3) )
				        map_mtxV(2,i) = map_mtxV(2,i) + (map_dtime/4.0) * &
				        	(-map_dc(1)*H_dot(i,1) + map_dc(3)*H_dot(i,3) - map_dc(2)*H_dot(i,4) )
				        map_mtxV(1,i) = map_mtxV(1,i) + (map_dtime/2.0) * &
				        	( map_dc(1)*H_dot(i,2) + map_dc(2)*H_dot(i,3) + map_dc(3)*H_dot(i,4) )	
				    end do
		    end select
		!-- IMAGINARY dynamics
		elseif(map_imag .eq. 1) then
		    select case (map_model)
		    	!-- 4-valuables model with normalization
			    case (1)
				    do i=1, map_mtxsize
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/2) * H_dot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/2) * H_dot(i,4)
					    map_mtxV(3,i) = map_mtxV(3,i) + (dtime)   * H_dot(i,1)
					    map_mtxV(4,i) = map_mtxV(4,i) + (dtime)   * H_dot(i,2)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime/2) * H_dot(i,3)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/2) * H_dot(i,4)
				    end do
				!-- 2-valuables model with normalization ( natural expression of Schrodinger Equation )
			    case (2)
				    do i=1,map_mtxsize
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/2) * H_dot(i,2)
					    map_mtxV(1,i) = map_mtxV(1,i) - (dtime)   * H_dot(i,1)
					    map_mtxV(2,i) = map_mtxV(2,i) - (dtime/2) * H_dot(i,2)
				    end do
				!-- 2-valuables model with parameter level off
				case (3)
					Elevel1 = dot_product( map_mtxV(1,:), dot_product(map_mtxH, map_mtxV(1,:)) )
			    	Elevel2 = dot_product( map_mtxV(2,:), dot_product(map_mtxH, map_mtxV(2,:)) )
			    	map_mtxV(2,:) = map_mtxV(2,:) - (dtime/2) * ( H_dot(i,2) - Elevel2 * map_mtxV(2,:) )
			    	map_mtxV(1,:) = map_mtxV(1,:) - (dtime) * ( H_dot(i,1) - Elevel1 * map_mtxV(1,:) )
			    	map_mtxV(2,:) = map_mtxV(2,:) - (dtime/2) * ( H_dot(i,2) - Elevel2 * map_mtxV(2,:) )
			    case default
				    print *, "not support now!"
				    stop
		    end select
		!-- evaluate the matrix of [ exp^(-H) ] 
		elseif(map_imag .eq. 2) then
			map_mtxO(:,:) = map_mtxO(:,:) - (dtime) * dot_product(map_mtxH, map_mtxO) 
		else
		    print *, "map_imag should be 0 or 1"
		    stop
		end if
	end subroutine map_propagator
	
	
	!-- initializer for time t=0
	!------ note that:
	!------ initialization calssified by key < map_initmd >
	subroutine map_initializer()
	implicit none
		integer :: i,j
		real    :: irand, temp1
		real(8), dimension(map_mtxsize) :: rtpop
		
		!-- differet model, different propagator
		!------ note that:
		!------ 1) map_mtxV(1,:)   stands for x
		!------ 1) map_mtxV(2,:)   stands for Px
		!------ 1) map_mtxV(3,:)   stands for y
		!------ 1) map_mtxV(4,:)   stands for Py
		
		!-- if needed
		call init_seed()
		
		!-- allocate root-population array
		if(map_istat .eq. 0) then ! population state =o means to be auto-arranged
		    temp1 = map_mtxsize*(2*map_mtxsize+1)*(map_mtxsize+1) / 6.0
		    temp1 = sqrt(temp1)
		    do i=1,map_mtxsize
		        rtpop(i) = real(i)/temp1
		    enddo
		else
		    rtpop(:) = 0
		    rtpop(map_istat) = 1
		endif
		
		select case (map_initmd)
		    case (1) ! map_model 1,3,4,5,6,7
		        do i=1,map_mtxsize
		            map_mtxV(1,i) = dcos(map_angle) * rtpop(i)
				    map_mtxV(2,i) = dsin(map_angle) * rtpop(i)
				    map_mtxV(3,i) = -dsin(map_angle) * rtpop(i)
				    map_mtxV(4,i) = dcos(map_angle) * rtpop(i)
			    end do
			case (2) ! map_model 2
			    do i=1,map_mtxsize
		            map_mtxV(1,i) = dsqrt(2._8) * dcos(map_angle) * rtpop(i)
				    map_mtxV(2,i) = dsqrt(2._8) * dsin(map_angle) * rtpop(i)
			    end do
			case (3) ! map_model 8, pseudo-rotation model
			    call random_number(irand)
			    map_cos = abs(cos(irand*map_2pi))
			    map_sin = abs(sin(irand*map_2pi))
			    map_cos = cos(map_pi*0.25)
			    map_sin = sin(map_pi*0.25)
			    print *, 'now my mapping angle (cos,sin) : ', map_cos, map_sin
			    do i=1,map_mtxsize
				    map_mtxV(1,i) = dcos(map_angle) * rtpop(i)
				    map_mtxV(2,i) = dsin(map_angle) * rtpop(i)
				    map_mtxV(3,i) = -dsin(map_angle) * rtpop(i)
				    map_mtxV(4,i) = dcos(map_angle) * rtpop(i)
			    end do
			case (4) ! map_model 9, spin-mapping model
			    do i=1,map_mtxsize
		            call random_number(irand)
		            if(i.eq. map_istat) then
		                !-- revise this part
					    map_mtxV(1,i) = irand*map_2pi
					    map_mtxV(2,i) = 0.5
				    else
					    map_mtxV(1,i) = irand*map_2pi
					    map_mtxV(2,i) = -0.5
				    end if
		        end do
		        print *, map_mtxV(1,:)
		    case (5) ! map_model 10, quaternion mapping
		        do i=1,map_mtxsize
		            map_mtxV(1,i) = dcos(map_angle)*dcos(map_angle) * rtpop(i)
				    map_mtxV(2,i) = dsin(map_angle)*dcos(map_angle) * rtpop(i)
				    map_mtxV(3,i) = dcos(map_angle)*dsin(map_angle) * rtpop(i)
				    map_mtxV(4,i) = dsin(map_angle)*dsin(map_angle) * rtpop(i)
				end do
		    case default
		        print *, "initmd wrong"
		        stop
		end select
	end subroutine map_initializer



	!-- analysis for population of the state
	!------ note that:
	!------ there are only tree kinds of analysis ways
	!------ 1) case 1 and case 6 are same
	!------ 2) case 2
	!------ 3) case 3, case 4, case 5 are same
	subroutine map_analyzer
	implicit none
	    real(8) :: pop, sqrtpop
		integer :: i
		
		!-- differet model, different propagator
		!------ note that:
		!------ 1) map_mtxV(1,:)   stands for x
		!------ 1) map_mtxV(2,:)   stands for Px
		!------ 1) map_mtxV(3,:)   stands for y
		!------ 1) map_mtxV(4,:)   stands for Py
		select case (map_analmd)
			case (1) ! model 1,6
				do i=1,map_mtxsize
					map_mtxN(i) = map_mtxV(1,i)*map_mtxV(4,i) - map_mtxV(3,i)*map_mtxV(2,i)				
				end do
			case (2) ! model 2
				do i=1,map_mtxsize
					map_mtxN(i) = 0.5_8 * ( map_mtxV(1,i)**2 + map_mtxV(2,i)**2 )	
				end do
			case (3) ! model 3,4,5,7
				do i=1,map_mtxsize
					map_mtxN(i) = 0.25_8 * ( (map_mtxV(1,i)+map_mtxV(4,i))**2 + (map_mtxV(3,i)-map_mtxV(2,i))**2 )	
				end do
		    case (4) ! model 8
			    do i=1,map_mtxsize
					map_mtxN(i) =  0.5_8 * ( (map_cos*map_mtxV(1,i) + map_sin*map_mtxV(4,i))**2 + &
					    (map_cos*map_mtxV(3,i) - map_sin*map_mtxV(2,i))**2 ) * 2*map_cos*map_sin &
					    + (0.5 - map_cos*map_sin) * (map_mtxV(1,i)*map_mtxV(4,i)-map_mtxV(2,i)*map_mtxV(3,i))				
				end do
		    case (5) ! model 9
				do i=1,map_mtxsize
					map_mtxN(i) =  (0.5 + map_mtxV(2,i) )
				end do
			case (6) ! model 10
			    do i=1,map_mtxsize
					map_mtxN(i) =  map_mtxV(1,i)**2 + map_mtxV(2,i)**2 + map_mtxV(3,i)**2 + map_mtxV(4,i)**2
				end do
			case default
			    print *, "analmd wrong"
			    stop
		end select
		
		!-- add scaler here
		if(map_imag.eq. 1) then
		    pop = sum( map_mtxN )
		    sqrtpop = sqrt(pop)
		    map_mtxN(:) = map_mtxN(:) / pop
		    map_mtxV(:,:) = map_mtxV(:,:) / sqrtpop
		endif
	end subroutine map_analyzer


	subroutine map_sampler()
	implicit none
		write(20,*) map_npass, map_mtxN
		write(21,*) map_npass, map_mtxV(1,1), map_mtxV(2,1)
		if(map_model.eq.2) then
    		write(22,*) map_npass, map_mtxV(1,2)/sqrt(2.0), map_mtxV(2,2)/sqrt(2.0)
        else
            write(22,*) map_npass, 0.5*(map_mtxV(1,2)+map_mtxV(4,2)), 0.5*(map_mtxV(2,2)-map_mtxV(3,2))
        end if
	end subroutine map_sampler

	!-- control the all propagation
	!------ note that; this binds the procedure "initializer, propagator, analyzer" all together
	subroutine map_controller
	implicit none
		integer :: i
		
		!-- sampling to pop.dat
		open(unit=20, file=trim('pop.dat'), status='replace')
		!-- coordinate data
		open(unit=21, file=trim('coo.dat'), status='replace')
		!-- project on 2nd state data
		open(unit=22, file=trim('pj2.dat'), status='replace')
		call map_initializer()
		
		map_npass = 0
		!-- REAL dynamics
		if(map_imag .eq. 0) then
			do i=1,map_nstep
				call map_propagator(map_dtime)
				call map_analyzer()
				map_npass = map_npass + 1
				if(mod(map_npass, map_nsamp) .eq. 0) then
					call map_sampler()
				end if
			enddo
		!-- IMAG dynamics
		elseif(map_imag .eq. 1) then
			!call map_eraser('el1.dat')
			!call map_eraser('el2.dat')
			!call map_eraser('el3.dat')
			!call map_renorm()
			do i=1,map_nstep
				call map_propagator(map_dtime)
				!call map_eraser('el1.dat')
				!call map_eraser('el2.dat')
				!call map_eraser('el3.dat')
				!call map_renorm()
				call map_analyzer()
				map_npass = map_npass + 1
				if(mod(map_npass, map_nsamp) .eq. 0) then
					call map_sampler()
				end if
			end do
		!-- Boltzmann OP
		else
			do i=1,map_nstep
				call map_propagator(map_beta/map_nstep)
			enddo
			open(unit=23, file=trim('bol.dat'), status='replace')
			
		endif
		
		close(unit=20)
		close(unit=21)
		close(unit=22)
		
		open(unit=20, file=trim('vs.dat'), status='replace')
		do i=1,map_mtxsize
			write(20,*) map_mtxV(:,i)
		end do
		close(unit=20)
		
	end subroutine map_controller

end module map4


!-- the main program of 
program main
use map4
implicit none
	call map_initmap()
	call map_controller()
	print *, "run down"
end program main

