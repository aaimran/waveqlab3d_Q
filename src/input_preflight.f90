module input_preflight

  use mpi
  use common, only : wp
  use datatypes, only : block_temp_parameters
  use diagnostics, only : diagnostic_list_t, DIAG_WARNING, DIAG_ERROR, &
       terminate_collectively
  use simulation_config, only : simulation_config_t
  use anelastic_q4_model, only : read_q4_parameters
  use anelastic_q8_model, only : read_q8_parameters
  use anelastic_fq8_model, only : read_fq8_parameters
  use decomposition_safety, only : stencil_requirements_t, get_stencil_requirements, &
       topology_fits, select_single_block_topology
  implicit none
  private

  public :: preflight_input

contains

  subroutine preflight_input(filename, config)
    character(*), intent(in) :: filename
    type(simulation_config_t), intent(out) :: config

    type(diagnostic_list_t) :: issues
    type(block_temp_parameters) :: btp(2)
    character(len=256) :: name, problem, response, plastic_model
    character(len=64) :: coupling, fd_type, mesh_source, type_of_mesh, material_source
    character(len=512) :: iomsg
    integer :: nblocks, nt, order, w_stride
    integer :: infile, stat, ierr, world_rank, world_size
    real(wp) :: CFL, t_final, topo
    logical :: w_fault, interpol, use_topography, mollify_source, valid

    namelist /problem_list/ name, problem, response, plastic_model, nblocks, &
         nt, CFL, coupling, fd_type, order, t_final, mesh_source, type_of_mesh, &
         material_source, interpol, w_stride, w_fault, use_topography, topo, &
         mollify_source
    namelist /block_list/ btp

    call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)

    if (world_rank == 0) then
       call set_problem_defaults(name, problem, response, plastic_model, nblocks, &
            nt, CFL, coupling, fd_type, order, t_final, mesh_source, type_of_mesh, &
            material_source, interpol, w_stride, w_fault, use_topography, topo, &
            mollify_source)

       open(newunit=infile, file=filename, status='old', action='read', &
            iostat=stat, iomsg=iomsg)
       if (stat /= 0) then
          call issues%add(DIAG_ERROR, 'CFG-FILE-001', &
               'Cannot open input file: '//trim(iomsg), section='input file', &
               field=trim(filename), suggestion='Check the path and read permissions.')
       else
          read(infile, nml=problem_list, iostat=stat, iomsg=iomsg)
          if (stat /= 0) then
             call issues%add(DIAG_ERROR, 'CFG-PROBLEM-001', &
                  'Cannot parse &problem_list: '//trim(iomsg), &
                  section='problem_list', suggestion='Fix namelist syntax or unknown fields.')
          else
             call validate_problem(response, nblocks, CFL, t_final, w_stride, &
                  fd_type, order, issues)
             if (trim(adjustl(response)) == 'anelastic-Q') then
                call issues%add(DIAG_WARNING, 'CFG-Q4-DEP-001', &
                     'Response anelastic-Q is deprecated and is normalized to anelastic-Q4.', &
                     section='problem_list', field='response', &
                     suggestion='Set response=''anelastic-Q4'' and use &anelastic_Q4_list.')
                response = 'anelastic-Q4'
             end if
             if (trim(adjustl(response)) == 'frequency-Q-8M') then
                call issues%add(DIAG_WARNING, 'CFG-FQ8-DEP-001', &
                     'Response frequency-Q-8M is deprecated and is normalized to anelastic-fQ8.', &
                     section='problem_list', field='response', &
                     suggestion='Set response=''anelastic-fQ8'' and use &anelastic_fQ8_list.')
                response = 'anelastic-fQ8'
             end if
          end if

          if (.not.issues%has_errors()) then
             rewind(infile)
             read(infile, nml=block_list, iostat=stat, iomsg=iomsg)
             if (stat /= 0) then
                call issues%add(DIAG_ERROR, 'CFG-BLOCK-001', &
                     'Cannot parse &block_list: '//trim(iomsg), &
                     section='block_list', suggestion='Define every active block and fix namelist syntax.')
             else
                call validate_blocks(btp, nblocks, fd_type, order, issues)
                if (nblocks == 2 .and. .not.issues%has_errors()) &
                     call validate_two_block_interface(btp, issues)
                if (.not.issues%has_errors()) &
                     call resolve_decomposition(btp, nblocks, world_size, fd_type, order, &
                          config, issues)
             end if
          end if

          if (.not.issues%has_errors() .and. trim(adjustl(response)) == 'anelastic-Q4') then
             call read_q4_parameters(infile, config%q4, stat, iomsg)
             if (stat /= 0) then
                call issues%add(DIAG_ERROR, 'CFG-Q4-001', trim(iomsg), &
                     section='anelastic_Q4_list', &
                     suggestion='Provide positive Qs0/Qp0 in &anelastic_Q4_list; legacy c is unsupported.')
             else
                config%has_q4 = .true.
             end if
          end if

          if (.not.issues%has_errors() .and. trim(adjustl(response)) == 'anelastic-Q8') then
             call read_q8_parameters(infile, config%q8, stat, iomsg)
             if (stat /= 0) then
                call issues%add(DIAG_ERROR, 'CFG-Q8-001', trim(iomsg), &
                     section='anelastic_Q8_list', &
                     suggestion='Provide positive Qs0/Qp0 and a supported coefficient setup.')
             else
                config%has_q8 = .true.
             end if
          end if
          if (.not.issues%has_errors() .and. trim(adjustl(response)) == 'anelastic-fQ8') then
             call read_fq8_parameters(infile, config%fq8, stat, iomsg)
             if (stat /= 0) then
                call issues%add(DIAG_ERROR, 'CFG-FQ8-001', trim(iomsg), &
                     section='anelastic_fQ8_list', &
                     suggestion='Use coarse_grain=2 with coefficient_method=''withers-2015'', '// &
                     'or use coarse_grain=0 with coefficient_method=''conventional-nnls''; '// &
                     'also provide Qs0/Qp0 >= 15 and valid gamma/frequencies.')
             else
                config%has_fq8 = .true.
             end if
          end if
          close(infile)
       end if

       valid = .not.issues%has_errors()
       if (issues%count() > 0) then
          if (valid) then
             call issues%print_summary('Input preflight passed with warnings: '//trim(filename))
          else
             call issues%print_summary('Input preflight failed: '//trim(filename))
          end if
       else
          write(*,'(A)') 'Input preflight passed: 0 errors, 0 warnings.'
       end if

       if (valid) then
          config%problem%name = name
          config%problem%problem = problem
          config%problem%response = trim(adjustl(response))
          config%problem%plastic_model = plastic_model
          config%problem%coupling = coupling
          config%problem%fd_type = fd_type
          config%problem%mesh_source = mesh_source
          config%problem%type_of_mesh = type_of_mesh
          config%problem%material_source = material_source
          config%problem%nblocks = nblocks
          config%problem%nt = nt
          config%problem%order = order
          config%problem%w_stride = w_stride
          config%problem%CFL = CFL
          config%problem%t_final = t_final
          config%problem%topo = topo
          config%problem%w_fault = w_fault
          config%problem%interpol = interpol
          config%problem%use_topography = use_topography
          config%problem%mollify_source = mollify_source
          allocate(config%blocks(nblocks))
          config%blocks = btp(1:nblocks)
          call print_decomposition(config)
       end if
    end if

    call MPI_Bcast(valid, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (.not.valid) call terminate_collectively(2)
    call broadcast_config(config, world_rank)
  end subroutine preflight_input


  subroutine print_decomposition(config)
    type(simulation_config_t), intent(in) :: config
    integer :: i

    write(*,'(A)') 'Resolved MPI decomposition:'
    do i = 1, config%problem%nblocks
       write(*,'(A,I0,A,3(I0,1X),A,3(I0,1X),A,I0,A,I0)') &
            '  block ', i, ': grid ', config%blocks(i)%nqrs, &
            ' topology ', config%process_dims(i,:), ' ranks ', &
            config%rank_begin(i), '-', config%rank_end(i)
    end do
    if (config%serial_shared_blocks) write(*,'(A)') &
         '  serial two-block mode: rank 0 owns both blocks'
  end subroutine print_decomposition


  subroutine set_problem_defaults(name, problem, response, plastic_model, nblocks, &
       nt, CFL, coupling, fd_type, order, t_final, mesh_source, type_of_mesh, &
       material_source, interpol, w_stride, w_fault, use_topography, topo, mollify_source)
    character(*), intent(out) :: name, problem, response, plastic_model
    character(*), intent(out) :: coupling, fd_type, mesh_source, type_of_mesh, material_source
    integer, intent(out) :: nblocks, nt, order, w_stride
    real(wp), intent(out) :: CFL, t_final, topo
    logical, intent(out) :: interpol, w_fault, use_topography, mollify_source

    name='default'; problem='TPV5'; response='elastic'; plastic_model='default'
    nblocks=2; nt=0; CFL=0.5_wp; coupling='locked'; fd_type='traditional'
    order=5; t_final=0.0_wp; mesh_source='compute'; type_of_mesh='cartesian'
    material_source='hardcode'; interpol=.false.; w_stride=1; w_fault=.true.
    use_topography=.false.; topo=1.0_wp; mollify_source=.false.
  end subroutine set_problem_defaults


  subroutine validate_problem(response, nblocks, CFL, t_final, w_stride, fd_type, order, issues)
    character(*), intent(in) :: response, fd_type
    integer, intent(in) :: nblocks, w_stride
    integer, intent(in) :: order
    real(wp), intent(in) :: CFL, t_final
    type(diagnostic_list_t), intent(inout) :: issues
    type(stencil_requirements_t) :: req

    select case (trim(adjustl(response)))
    case ('elastic','plastic','anelastic','low-pass','anelastic-Q','anelastic-Q4','anelastic-Q8', &
          'anelastic-Qf','anelastic-fQ8','constant-Q-4M','constant-Q-8M','frequency-Q-4M','frequency-Q-8M')
    case default
       call issues%add(DIAG_ERROR, 'CFG-PROBLEM-002', &
            'Unsupported response: '//trim(response), section='problem_list', field='response')
    end select
    if (nblocks < 1 .or. nblocks > 2) call issues%add(DIAG_ERROR, 'CFG-PROBLEM-003', &
         'nblocks must be 1 or 2.', section='problem_list', field='nblocks')
    if (CFL <= 0.0_wp) call issues%add(DIAG_ERROR, 'CFG-PROBLEM-004', &
         'CFL must be positive.', section='problem_list', field='CFL')
    if (t_final < 0.0_wp) call issues%add(DIAG_ERROR, 'CFG-PROBLEM-005', &
         't_final cannot be negative.', section='problem_list', field='t_final')
    if (w_stride < 1) call issues%add(DIAG_ERROR, 'CFG-PROBLEM-006', &
         'w_stride must be at least 1.', section='problem_list', field='w_stride')
    if (t_final == 0.0_wp) call issues%add(DIAG_WARNING, 'CFG-PROBLEM-007', &
         't_final is zero; initialization will run but no time steps are expected.', &
         section='problem_list', field='t_final')
    req = get_stencil_requirements(fd_type, order)
    if (.not.req%supported) call issues%add(DIAG_ERROR, 'CFG-FD-001', &
         'Unsupported finite-difference type/order combination: '//trim(fd_type), &
         section='problem_list', field='fd_type/order', &
         suggestion='Use traditional, upwind order 2-9, or a supported upwind_drp order.')
  end subroutine validate_problem


  subroutine validate_blocks(btp, nblocks, fd_type, order, issues)
    type(block_temp_parameters), intent(in) :: btp(2)
    integer, intent(in) :: nblocks
    character(*), intent(in) :: fd_type
    integer, intent(in) :: order
    type(diagnostic_list_t), intent(inout) :: issues
    integer :: i
    type(stencil_requirements_t) :: req

    req = get_stencil_requirements(fd_type, order)

    do i = 1, nblocks
       if (any(btp(i)%nqrs < 2)) call issues%add(DIAG_ERROR, 'CFG-BLOCK-002', &
            'Every active grid dimension must contain at least 2 points.', &
            section='block_list', field='nqrs', block_id=i)
       if (any(btp(i)%bqrs <= btp(i)%aqrs)) call issues%add(DIAG_ERROR, 'CFG-BLOCK-003', &
            'Every bqrs coordinate must exceed the corresponding aqrs coordinate.', &
            section='block_list', field='aqrs/bqrs', block_id=i)
       if (any(btp(i)%rho_s_p <= 0.0_wp)) call issues%add(DIAG_ERROR, 'CFG-BLOCK-004', &
            'Density, Vs, and Vp must be positive.', section='block_list', &
            field='rho_s_p', block_id=i)
       if (btp(i)%npml < 0) call issues%add(DIAG_ERROR, 'CFG-BLOCK-005', &
            'npml cannot be negative.', section='block_list', field='npml', block_id=i)
       if (req%supported .and. any(btp(i)%nqrs < req%minimum_global_points)) &
            call issues%add(DIAG_ERROR, 'CFG-STENCIL-001', &
            'Grid is too small for nonoverlapping physical-boundary closures.', &
            section='block_list', field='nqrs', block_id=i, &
            suggestion='Increase every grid dimension to the operator global minimum.')
       if (btp(i)%npml > 0) then
          if (any(btp(i)%nqrs <= btp(i)%npml*(merge(1,0,btp(i)%pml_lqrs) + &
                                              merge(1,0,btp(i)%pml_rqrs)))) &
               call issues%add(DIAG_ERROR, 'CFG-PML-001', &
               'PML layers leave no non-PML core in at least one direction.', &
               section='block_list', field='npml/pml_lqrs/pml_rqrs', block_id=i)
       end if
    end do
  end subroutine validate_blocks


  subroutine validate_two_block_interface(btp, issues)
    type(block_temp_parameters), intent(in) :: btp(2)
    type(diagnostic_list_t), intent(inout) :: issues
    real(wp) :: scale, tolerance

    scale = max(1.0_wp, maxval(abs([btp(1)%aqrs, btp(1)%bqrs, &
                                    btp(2)%aqrs, btp(2)%bqrs])))
    tolerance = 1000.0_wp*epsilon(1.0_wp)*scale

    if (any(btp(1)%nqrs(2:3) /= btp(2)%nqrs(2:3))) then
       call issues%add(DIAG_ERROR, 'CFG-IFACE-001', &
            'The two blocks must have identical r/s (Y-Z) grid counts.', &
            section='block_list', field='btp%nqrs(2:3)', &
            suggestion='Keep Y and Z counts equal; only the q/X count may differ.')
    end if
    if (any(abs(btp(1)%aqrs(2:3)-btp(2)%aqrs(2:3)) > tolerance) .or. &
        any(abs(btp(1)%bqrs(2:3)-btp(2)%bqrs(2:3)) > tolerance)) then
       call issues%add(DIAG_ERROR, 'CFG-IFACE-002', &
            'The two blocks must have identical r/s (Y-Z) physical extents.', &
            section='block_list', field='btp%aqrs/bqrs(2:3)')
    end if
    if (abs(btp(1)%bqrs(1)-btp(2)%aqrs(1)) > tolerance) then
       call issues%add(DIAG_ERROR, 'CFG-IFACE-003', &
            'Block 1 right-q coordinate must equal block 2 left-q coordinate.', &
            section='block_list', field='q interface coordinate')
    end if
    if (btp(1)%rqrs(1) /= 0 .or. btp(2)%lqrs(1) /= 0) then
       call issues%add(DIAG_ERROR, 'CFG-IFACE-004', &
            'Internal q faces must both use interface boundary code 0.', &
            section='block_list', field='btp(1)%rqrs(1), btp(2)%lqrs(1)')
    end if
    if (btp(1)%pml_rqrs(1) .or. btp(2)%pml_lqrs(1)) then
       call issues%add(DIAG_ERROR, 'CFG-IFACE-005', &
            'PML cannot be enabled on either side of the internal interface.', &
            section='block_list', field='internal q-face PML')
    end if
  end subroutine validate_two_block_interface


  subroutine resolve_decomposition(btp, nblocks, world_size, fd_type, order, config, issues)
    type(block_temp_parameters), intent(in) :: btp(2)
    integer, intent(in) :: nblocks, world_size
    character(*), intent(in) :: fd_type
    integer, intent(in) :: order
    type(simulation_config_t), intent(inout) :: config
    type(diagnostic_list_t), intent(inout) :: issues
    integer :: dims(3), pq1, pq2, pr, ps, tangential, qsum, minimum_owned
    integer :: best_dims(2,3), best_sizes(2)
    real(wp) :: work(2), load1, load2, score, best_score
    type(stencil_requirements_t) :: req
    logical :: found

    req = get_stencil_requirements(fd_type, order)
    minimum_owned = req%minimum_owned_points

    config%process_dims = 1
    config%block_sizes = 0
    config%rank_begin = 0
    config%rank_end = -1
    config%serial_shared_blocks = .false.

    if (nblocks == 1) then
       call select_single_block_topology(btp(1)%nqrs, world_size, minimum_owned, dims, found)
       if (.not.found) then
          call issues%add(DIAG_ERROR, 'CFG-DECOMP-002', &
               'No Cartesian factorization satisfies the owned-subdomain bound.', &
               section='decomposition', block_id=1, &
               suggestion='Use fewer ranks or enlarge the grid.')
          return
       end if
       config%process_dims(1,:) = dims
       config%block_sizes(1) = world_size
       config%rank_begin(1) = 0
       config%rank_end(1) = world_size-1
       return
    end if

    if (world_size == 1) then
       config%process_dims(1,:) = 1
       config%process_dims(2,:) = 1
       config%block_sizes = 1
       config%rank_begin = 0
       config%rank_end = 0
       config%serial_shared_blocks = .true.
       return
    end if

    work(1) = real(product(btp(1)%nqrs), wp)
    work(2) = real(product(btp(2)%nqrs), wp)
    best_score = huge(1.0_wp)
    best_dims = 0
    best_sizes = 0

    do pr = 1, world_size
       do ps = 1, world_size/pr
          tangential = pr*ps
          if (mod(world_size, tangential) /= 0) cycle
          qsum = world_size/tangential
          if (qsum < 2) cycle
          do pq1 = 1, qsum-1
             pq2 = qsum-pq1
             dims = [pq1, pr, ps]
             if (.not.topology_fits(btp(1)%nqrs, dims, minimum_owned)) cycle
             dims = [pq2, pr, ps]
             if (.not.topology_fits(btp(2)%nqrs, dims, minimum_owned)) cycle
             load1 = work(1)/real(pq1*tangential,wp)
             load2 = work(2)/real(pq2*tangential,wp)
             score = abs(load1-load2)/max(load1,load2)
             if (score < best_score) then
                best_score = score
                best_dims(1,:) = [pq1,pr,ps]
                best_dims(2,:) = [pq2,pr,ps]
                best_sizes = [pq1*tangential,pq2*tangential]
             end if
          end do
       end do
    end do

    if (best_sizes(1) == 0) then
       call issues%add(DIAG_ERROR, 'CFG-DECOMP-101', &
            'No asymmetric two-block topology satisfies the shared Y-Z topology and 20-point bound.', &
            section='decomposition', &
            suggestion='Use fewer ranks, enlarge a block, or choose a factorable world size.')
       return
    end if

    config%process_dims = best_dims
    config%block_sizes = best_sizes
    config%rank_begin = [0,best_sizes(1)]
    config%rank_end = [best_sizes(1)-1,world_size-1]
  end subroutine resolve_decomposition


  subroutine broadcast_config(config, world_rank)
    type(simulation_config_t), intent(inout) :: config
    integer, intent(in) :: world_rank
    integer :: ierr, i, nblocks

    nblocks = config%problem%nblocks
    call MPI_Bcast(nblocks, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    config%problem%nblocks = nblocks
    if (world_rank /= 0) allocate(config%blocks(nblocks))

    call bcast_chars(config%problem%name)
    call bcast_chars(config%problem%problem)
    call bcast_chars(config%problem%response)
    call bcast_chars(config%problem%plastic_model)
    call bcast_chars(config%problem%coupling)
    call bcast_chars(config%problem%fd_type)
    call bcast_chars(config%problem%mesh_source)
    call bcast_chars(config%problem%type_of_mesh)
    call bcast_chars(config%problem%material_source)
    call MPI_Bcast(config%problem%nt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%problem%order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%problem%w_stride, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%problem%CFL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%problem%t_final, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%problem%topo, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%problem%w_fault, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%problem%interpol, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%problem%use_topography, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%problem%mollify_source, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%process_dims, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%block_sizes, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%rank_begin, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%rank_end, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%serial_shared_blocks, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    do i = 1, nblocks
       call broadcast_block(config%blocks(i))
    end do
    call MPI_Bcast(config%has_q4, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (config%has_q4) call broadcast_q4(config)
    call MPI_Bcast(config%has_q8, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (config%has_q8) call broadcast_q8(config)
    call MPI_Bcast(config%has_fq8, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (config%has_fq8) call broadcast_fq8(config)
  end subroutine broadcast_config


  subroutine broadcast_block(block)
    type(block_temp_parameters), intent(inout) :: block
    integer :: ierr
    call MPI_Bcast(block%nqrs, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%aqrs, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%bqrs, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%lc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%rc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%rho_s_p, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%mu_beta_eta, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%lqrs, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%rqrs, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call bcast_chars(block%profile_type); call bcast_chars(block%profile_path)
    call bcast_chars(block%material_path(1)); call bcast_chars(block%material_path(2))
    call bcast_chars(block%material_path(3)); call bcast_chars(block%topography_type)
    call bcast_chars(block%topography_path)
    call MPI_Bcast(block%pml_lqrs, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%pml_rqrs, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%npml, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(block%faultsize, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end subroutine broadcast_block


  subroutine broadcast_q4(config)
    type(simulation_config_t), intent(inout) :: config
    integer :: ierr
    call bcast_chars(config%q4%weight_method)
    call MPI_Bcast(config%q4%Qs0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%q4%Qp0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%q4%fref, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%q4%fmin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%q4%fmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine broadcast_q4


  subroutine broadcast_q8(config)
    type(simulation_config_t), intent(inout) :: config
    integer :: ierr
    call bcast_chars(config%q8%weight_method)
    call MPI_Bcast(config%q8%Qs0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%q8%Qp0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%q8%fref, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%q8%fmin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%q8%fmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine broadcast_q8

  subroutine broadcast_fq8(config)
    type(simulation_config_t), intent(inout) :: config
    integer :: ierr
    call bcast_chars(config%fq8%coefficient_method)
    call bcast_chars(config%fq8%weight_policy)
    call MPI_Bcast(config%fq8%coarse_grain, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%fq8%Qs0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%fq8%Qp0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%fq8%gamma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%fq8%f_transition, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(config%fq8%fref, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine broadcast_fq8


  subroutine bcast_chars(value)
    character(*), intent(inout) :: value
    integer :: ierr
    call MPI_Bcast(value, len(value), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  end subroutine bcast_chars

end module input_preflight
