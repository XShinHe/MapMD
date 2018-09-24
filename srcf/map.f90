!*
!--- Copyright by --- XShinHe <1500011805@pku.edu.cn>
!------- Date 2018. 08
!--- Acknowledgement to Liu group, of PKU
!*

module map
implicit none
	!-- dimension of the matrix
	integer :: map_mtxsize
	!-- allocatable object of H-matrix
	real(8), dimension(:,:), allocatable :: map_mtxH ! Hamiltonian
	real(8), dimension(:,:), allocatable :: map_mtxV ! kinds of valuables
	real(8), dimension(:), allocatable :: map_mtxN   ! Populations
	
	!-- some arguments for different models
	integer :: map_model
	real(8) :: map_angle
	integer :: map_vardim
	
	!-- modification of Hamiltonian
	real(8) :: map_Hnondiag = -1.0
	integer :: map_arrange = 0
	
	!-- general for (real-time) dynamcs simulation
	real(8) :: map_dtime
	integer :: map_nstep
	integer :: map_npass
	integer :: map_nsamp
	integer :: map_istat ! the initial populated state
	
	real(8) :: map_pi  = 3.14159265358979323846
	real(8) :: map_2pi = 2.0_8 * 3.14159265358979323846
	real(8) :: map_cos
	real(8) :: map_sin
	
	!-- general for (imaginary) dynamics simulations & and thermstat-arguments
	real(8) :: map_temp
	real(8) :: map_beta
	! real(8), dimension(4,2) :: langc
	!-------- :: UNUSUAL THERMOSTAT
	

contains
	!-- defination of dot operator of non-diagonal elements
	!------ note that: matrix H and vector x should have same size in one direction, here we don't verify anymore
	!------ note that: H_dot = H_ndot + H_ddot     ( nondiagonal + diagonal )
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
	
	real(8) function H_xdot(n,flag) ! half-diagonal
		integer, intent(in) :: n, flag
		H_xdot = dot_product( map_mtxH(n,:), map_mtxV(flag,:) )
		H_xdot = H_xdot - 0.5_8 * map_mtxH(n,n)*map_mtxV(flag,n)
	end function H_xdot
	
	real(8) function spinm(n)
	    integer, intent(in) :: n
	    real(8) :: S2
	    S2 = 0.75
	    spinm = dsqrt( S2 - map_mtxV(2,n)**2 )
	end function spinm
	
	real(8) function cosd(i,j)
	    integer, intent(in) :: i,j
	    cosd = dcos( map_mtxV(2,i) - map_mtxV(1,j) )
	end function cosd
	
	real(8) function sind(i,j)
	    integer, intent(in) :: i,j
	    sind = dsin( map_mtxV(1,i) - map_mtxV(1,j) )
	end function sind
	
	
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

	!-- read from maatrx-file and initialize the Hamiltonian matirx
	!------ the default read-file is 'Hsave.dat'
	subroutine map_initmap()
	implicit none
		character(len=20), dimension(2) :: pairs
		integer :: i,j
		integer :: my_count, my_iostat
		logical :: my_exist
		
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
		                if(map_model .eq. 2) then
		                	map_vardim = 2
		                else
		                	map_vardim = 4
		                end if
		            case ('phsangle')
		            	read(pairs(2),*) map_angle
		            case ('mapdtime')
		            	read(pairs(2),*) map_dtime
		            case ('mapbeta')
		            	read(pairs(2),*) map_beta
		            case ('mapnstep')
		            	read(pairs(2),*) map_nstep
		            case ('mapnsamp')
		            	read(pairs(2),*) map_nsamp
		            case ('inistate')
		            	read(pairs(2),*) map_istat
		            case ('hnondiag')
		            	read(pairs(2),*) map_Hnondiag
		            case ('maparrange')
		            	read(pairs(2),*) map_arrange
		            case default
		            	print *, "warning, arguments in 'map.rc' doesn't match!"
		            	my_count = my_count - 1
		        end select
		    end do
		else
			print *, "error, 'map.rc' doesn't exist!"
			stop
		end if
		
		!-- verify the completeness of arguments
		if(my_count.eq. 9) then
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
	end subroutine map_initmap
	
	!-- propagate the assistant valuables
	!------ note MAYBEERROR: the discription of 'step1, step2, step3' are binded in a loop, or seperate to different loops ?
	subroutine map_propagator(dtime)
	implicit none
		real(8), intent(in) :: dtime
		integer :: i,j
		real    :: irand, x, y, c, s, cc, ss, sc
		
		!-- differet model, different propagator
		!------ note that:
		!------ 1) map_mtxV(1,:)   stands for x
		!------ 1) map_mtxV(2,:)   stands for Px
		!------ 1) map_mtxV(3,:)   stands for y
		!------ 1) map_mtxV(4,:)   stands for Py
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
			    irand = 20
				do i=1,map_mtxsize
				    !call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*H_ndot(i,3) - y*H_ndot(i,2) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * (-x*H_ddot(i,3) + y*H_ddot(i,2) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*H_ndot(i,3) - y*H_ndot(i,2) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
				!end do
				!do i=1,map_mtxsize
				    !call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * ( x*H_ddot(i,3) - y*H_ddot(i,2) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * ( y*H_ndot(i,2) - x*H_ndot(i,3) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(3,i) = map_mtxV(3,i) + (dtime/4) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(3,i) = map_mtxV(3,i) + (dtime/2) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(1,i) = map_mtxV(1,i) + (dtime/4) * ( x*H_ndot(i,2) - y*H_ndot(i,3) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(3,i) = map_mtxV(3,i) + (dtime/4) * ( y*H_ndot(i,4) + x*H_ndot(i,1) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(1,i) = map_mtxV(1,i) - (dtime/4) * ( x*H_ddot(i,3) - y*H_ddot(i,2) )
				!end do
				!do i=1,map_mtxsize
				    !call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*H_ndot(i,3) - y*H_ndot(i,2) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(4,i) = map_mtxV(4,i) + (dtime/4) * (-x*H_ddot(i,3) + y*H_ddot(i,2) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ndot(i,1) + y*H_ndot(i,4) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(4,i) = map_mtxV(4,i) - (dtime/8) * ( x*H_ndot(i,3) - y*H_ndot(i,2) )
					!call random_number(irand)
				    x = 2*irand
				    y = 2 - x
					map_mtxV(2,i) = map_mtxV(2,i) - (dtime/8) * ( x*H_ddot(i,1) + y*H_ddot(i,4) )
				end do
		    case (8)
		        c = sqrt(2.0)*map_cos
		        s = sqrt(2.0)*map_sin
		        cc = sqrt(2.0)*map_cos*map_sin
		        ss = sqrt(2.0)*map_cos*map_sin
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
			!-- map_mtxV(1,i):    m
			!-- map_mtxV(2,i):    q
			case (9)
				do i=1, map_mtxsize
				    map_mtxV(1,i) = map_mtxV(1,i) + (dtime) * (1.0) * map_mtxH(i,i)
				    do j=1, map_mtxsize
				        !-- cancel the next line
				        if(i.eq.j) cycle
				        map_mtxV(1,i) = map_mtxV(1,i) - (dtime) * 2* (map_mtxV(2,i)/spinm(i)) * spinm(j) * cosd(i,j) *map_mtxH(i,j)
				    end do
				end do
				do i=1, map_mtxsize
				    do j=1, map_mtxsize
				        if(i.eq.j) cycle
				        map_mtxV(2,i) = map_mtxV(2,i) + (dtime) * 2*spinm(i) * spinm(j) * sind(i,j) *map_mtxH(i,j)
				    end do
				end do
		end select
	end subroutine map_propagator
	
	
	!-- initializer foor t=0
	!------ note that:
	!------ there are actually two kinds of initialization
	!------ 1) two-dims valuables  : see case 2
	!------ 2) four-dims valuables : see case 1,3,4,5,6
	subroutine map_initializer(ketn)
	implicit none
		!-- initialize on the n-th state
		integer, intent(in) :: ketn
		integer :: i
		real    :: irand, temp1
		real(8), dimension(map_mtxsize) :: pop
		
		!-- differet model, different propagator
		!------ note that:
		!------ 1) map_mtxV(1,:)   stands for x
		!------ 1) map_mtxV(2,:)   stands for Px
		!------ 1) map_mtxV(3,:)   stands for y
		!------ 1) map_mtxV(4,:)   stands for Py
		
		!-- if needed
		call init_seed()
		
		if(map_arrange.eq. 1) then
		    temp1 = map_mtxsize*(2*map_mtxsize+1)*(map_mtxsize+1) / 6.0
		    temp1 = sqrt(temp1)
		    do i=1,map_mtxsize
		        pop(i) = real(i)/temp1
		    enddo
		else
		    pop(:) = 0
		    pop(ketn) = 1
		endif
		
		print *, pop
		
		!-- for case(map_model) 2
		if(map_model .eq. 2) then
		    do i=1,map_mtxsize
		        map_mtxV(1,i) = dsqrt(2._8) * dcos(map_angle) * pop(i)
				map_mtxV(2,i) = dsqrt(2._8) * dsin(map_angle) * pop(i)
			end do
	    !-- for case map_model 8
	    else if (map_model .eq. 8) then
	        call random_number(irand)
			map_cos = abs(cos(irand*map_2pi))
			map_sin = abs(sin(irand*map_2pi))
			map_cos = cos(map_pi*0.25)
			map_sin = sin(map_pi*0.25)
			print *, 'now my mapping angle (cos,sin) : ', map_cos, map_sin
			do i=1,map_mtxsize
				map_mtxV(1,i) = dcos(map_angle) * pop(i)
				map_mtxV(2,i) = dsin(map_angle) * pop(i)
				map_mtxV(3,i) = -dsin(map_angle) * pop(i)
				map_mtxV(4,i) = dcos(map_angle) * pop(i)
			end do
	    else if (map_model .eq. 9) then
		    do i=1,map_mtxsize
		        call random_number(irand)
		        if(i.eq.ketn) then
		            !-- revise this part
					map_mtxV(1,i) = irand*map_2pi
					map_mtxV(2,i) = 0.5
				else
					map_mtxV(1,i) = irand*map_2pi
					map_mtxV(2,i) = -0.5
				end if
		    end do
		    print *, map_mtxV(1,:)
		!-- for case(map_model) 1,3,4,5,6,7
		else
		    do i=1,map_mtxsize
		        map_mtxV(1,i) = dcos(map_angle) * pop(i)
				map_mtxV(2,i) = dsin(map_angle) * pop(i)
				map_mtxV(3,i) = -dsin(map_angle) * pop(i)
				map_mtxV(4,i) = dcos(map_angle) * pop(i)
			end do
		end if
	end subroutine map_initializer



	!-- analysis for population of the state
	!------ note that:
	!------ there are only tree kinds of analysis ways
	!------ 1) case 1 and case 6 are same
	!------ 2) case 2
	!------ 3) case 3, case 4, case 5 are same
	subroutine map_analyzer
	implicit none
		integer :: i
		
		!-- differet model, different propagator
		!------ note that:
		!------ 1) map_mtxV(1,:)   stands for x
		!------ 1) map_mtxV(2,:)   stands for Px
		!------ 1) map_mtxV(3,:)   stands for y
		!------ 1) map_mtxV(4,:)   stands for Py
		select case (map_model)
			case (1)
				do i=1,map_mtxsize
					map_mtxN(i) = map_mtxV(1,i)*map_mtxV(4,i) - map_mtxV(3,i)*map_mtxV(2,i)				
				end do
			case (2)
				do i=1,map_mtxsize
					map_mtxN(i) = 0.5_8 * ( map_mtxV(1,i)**2 + map_mtxV(2,i)**2 )	
				end do
			case (3)
				do i=1,map_mtxsize
					map_mtxN(i) = 0.25_8 * ( (map_mtxV(1,i)+map_mtxV(4,i))**2 + (map_mtxV(3,i)-map_mtxV(2,i))**2 )	
				end do
			case (4)
				do i=1,map_mtxsize
					map_mtxN(i) = 0.25_8 * ( (map_mtxV(1,i)+map_mtxV(4,i))**2 + (map_mtxV(3,i)-map_mtxV(2,i))**2 )	
				end do
			case (5)
				do i=1,map_mtxsize
					map_mtxN(i) = 0.25_8 * ( (map_mtxV(1,i)+map_mtxV(4,i))**2 + (map_mtxV(3,i)-map_mtxV(2,i))**2 )	
				end do
			case (6)
				do i=1,map_mtxsize
					map_mtxN(i) = map_mtxV(1,i)*map_mtxV(4,i) - map_mtxV(3,i)*map_mtxV(2,i)				
				end do
			case (7)
			    do i=1,map_mtxsize
					map_mtxN(i) =  0.25_8 * ( (map_mtxV(1,i)+map_mtxV(4,i))**2 + (map_mtxV(3,i)-map_mtxV(2,i))**2 )				
				end do
		    case (8)
			    do i=1,map_mtxsize
					map_mtxN(i) =  0.5_8 * ( (map_cos*map_mtxV(1,i) + map_sin*map_mtxV(4,i))**2 + &
					    (map_cos*map_mtxV(3,i) - map_sin*map_mtxV(2,i))**2 ) * 2*map_cos*map_sin &
					    + (0.5 - map_cos*map_sin) * (map_mtxV(1,i)*map_mtxV(4,i)-map_mtxV(2,i)*map_mtxV(3,i))				
				end do
		    case (9)
				do i=1,map_mtxsize
					map_mtxN(i) =  (0.5 + map_mtxV(2,i) )
				end do
		end select
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
	!------ note that; this binds the procedure "initializer, propagator, annlyzer" all together
	subroutine map_controller
	implicit none
		integer :: i
		
		!-- sampling to pop.dat
		open(unit=20, file=trim('pop.dat'), status='replace')
		open(unit=21, file=trim('coo.dat'), status='replace')
		open(unit=22, file=trim('pj2.dat'), status='replace')
		call map_initializer(map_istat)
		map_npass = 0
		do i=1,map_nstep
			call map_propagator(map_dtime)
			!call map_thermostat()
			call map_analyzer()
			map_npass = map_npass + 1
			if(mod(map_npass, map_nsamp) .eq. 0) then
				call map_sampler()
			end if
		end do
		close(unit=20)
		close(unit=21)
		close(unit=22)
		
		open(unit=20, file=trim('vs.dat'), status='replace')
		do i=1,map_mtxsize
			write(20,*) map_mtxV(:,i)
		end do
		close(unit=20)
		
	end subroutine map_controller


    real(8) function map_pop(i)
	implicit none
		integer, intent(in) :: i
		
		!-- differet model, different propagator
		!------ note that:
		!------ 1) map_mtxV(1,:)   stands for x
		!------ 1) map_mtxV(2,:)   stands for Px
		!------ 1) map_mtxV(3,:)   stands for y
		!------ 1) map_mtxV(4,:)   stands for Py
		if(map_model.eq. 1 .or. map_pop .eq. 6) then
		    map_pop = map_mtxV(1,i)*map_mtxV(4,i) - map_mtxV(3,i)*map_mtxV(2,i)
		else if (map_model .eq. 2) then
		    map_pop = 0.5_8 * ( map_mtxV(1,i)**2 + map_mtxV(2,i)**2 )
		else if (map_model < 8) then
		    map_pop = 0.25_8 * ( (map_mtxV(1,i)+map_mtxV(4,i))**2 + (map_mtxV(3,i)-map_mtxV(2,i))**2 )
		else
		    print *, 'error'
		    stop
		end if
	end function map_pop
	
	
    !-- thermostat Anderson
    subroutine map_thermostat()
    implicit none
        integer :: i,j
        real(8) :: irand, pi, pj, pij, deltah
        
        if(map_model.ne.2) then
            print *, 'thermostat only for model 2'
            stop
        end if
        do i=1,map_mtxsize
            do j=i+1,map_mtxsize
                call random_number(irand)
                if( irand > 0.9 ) then
                    deltah = map_beta*(map_mtxH(i,i) - map_mtxH(j,j))
                    pij = map_pop(i) + map_pop(j)
                    pi = 1.0/(1.0+exp(deltah))
                    pj = 1 - pi
                    call random_number(irand)
                    map_mtxV(1,i) = sqrt(2.0*pi*pij) * cos(irand*map_2pi)
                    map_mtxV(2,i) = sqrt(2.0*pi*pij) * sin(irand*map_2pi)
                    call random_number(irand)
                    map_mtxV(1,j) = sqrt(2.0*pj*pij) * cos(irand*map_2pi)
                    map_mtxV(2,j) = sqrt(2.0*pj*pij) * sin(irand*map_2pi)
                end if
            end do
        end do
    end subroutine map_thermostat

end module map


!-- the main program of 
program main
use map
implicit none
	call map_initmap()
	call map_controller()
	print *, "run down"
end program main


