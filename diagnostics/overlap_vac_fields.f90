program overlap_vac_fields

  use mod_parameters
  use data_structure
  use phys_module
  use basis_at_gaussian
  use mpi_mod
  use mod_import_restart
  use mod_export_restart
  use mod_boundary
  use mod_interp
  use mod_poiss
  use tr_module

  implicit none

  type (type_node_list)   ,     pointer :: node_list
  type (type_node_list)   ,     pointer :: aux_node_list
  type (type_element_list),     pointer :: element_list
  type (type_bnd_element_list), pointer :: bnd_elm_list    
  type (type_bnd_node_list),    pointer :: bnd_node_list 

  integer   :: ierr, my_id, i_harm, iter, i
  integer   :: n_cpu, n_cpu_n, n_cpu_master, m_cpu, n_masters, n_cpu_trans, my_id_trans
  integer   :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS
  integer   :: required,provided,StatInfo, resultlength
  integer*4 :: rank, comm_size 
  real*8    :: diff
  character(len=MPI_MAX_PROCESSOR_NAME) :: name

  allocate(node_list)
  allocate(aux_node_list) 
  allocate(element_list)
  allocate(bnd_elm_list)
  allocate(bnd_node_list)

  ! Initialization of memory, parallel threads etc
#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif

  call MPI_Init_thread(required, provided, StatInfo)

  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
  
  ! --- Determine number of MPI procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  n_cpu = comm_size
  
  ! --- Determine ID of each MPI proc
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  my_id = rank
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  CALL MPI_GET_PROCESSOR_NAME (name,resultlength,ierr)
  write(*,'(A,I5,2A)') '  #MPI id, ProcessorName ', rank, ': ', name
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call tr_meminit(my_id,n_cpu)

  ! Read JOREK input file
  call initialise_parameters(my_id, "__NO_FILENAME__")

  ! Initialize basis functions
  call initialise_basis    

  ! Initialize mode and mode_type arrays
  call det_modes()
  
  ! Read JOREK restart file
  call import_restart(node_list,  element_list, 'jorek_restart', rst_format, ierr, .true.)

  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

  ! Solve Poisson equation with boundary conditions
  i_harm = 2

  ! pure vac fields (no plasma current) use type option -3
  call poisson(my_id,-3,node_list,element_list,bnd_node_list,bnd_elm_list,var_psi,var_psi,i_harm, &
               0.d0,1.d0,.false.,1,-2.d0,.false.,.false.,1) ! These last parameters play no role here...   

  ! Export restart
  call export_restart(node_list, element_list, 'test')


end program overlap_vac_fields
