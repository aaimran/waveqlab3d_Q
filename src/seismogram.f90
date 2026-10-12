module seismogram

  use receiver_output_config, only: receiver_config_t, valid_station_name, lowercase
  use mpi
  use common, only : wp
  use datatypes, only : block_grid_t, seismogram_type
   use mpi3dbasic, only : is_master, boxed_lines
  implicit none

contains


  subroutine init_seismogram(input,S,name,G,cfg,nblocks)

    use diagnostics, only : fatal_local

    implicit none

    type(receiver_config_t),intent(in) :: cfg
    integer,intent(in) :: nblocks
    integer :: winner, candidate, ierr, bnd
    real(wp) :: best, local_distance
    integer,intent(in) :: input
    type(block_grid_t), intent(in) :: G
    type(seismogram_type),intent(inout) :: S
    character(*),intent(in) :: name

      integer :: mx,my,mz,px,py,pz
      logical :: output_exact_moment

      logical :: output_seismograms,output_fields_block1,output_fields_block2
      logical :: output_station_info
      logical :: output_station_mapping
      logical :: station_xyz_index
      logical :: station_number_in_list, station_number_in_filename
      logical :: station_use_block_subdirectories
      logical :: station_add_header, station_add_metadata
      logical :: any_station_outputs

      integer :: stride_fields,n,stat
      integer :: external_file_unit

      character(3) :: field
      character(32) :: xs, ys, zs
      character(32) :: block_str
      character(256) :: temp,filename
      character(256) :: station_list, station_list_file
      character(256) :: station_file_directory
      character(256) :: station_output_order
      character(256) :: common_stations_blocks
      character(256) :: station_info_lines(15)
      character(512) :: filepath
      character(512) :: block_output_directory
      character(512) :: station_output_directory

    !real(kind = wp) :: xmin, xmax, ymin, ymax, zmin, zmax



    mx = G%C%mq
    my = G%C%mr
    mz = G%C%ms
    px = G%C%pq
    py = G%C%pr
    pz = G%C%ps

    output_exact_moment=cfg%output_exact_moment
    output_seismograms=cfg%output_seismograms
    output_station_info=cfg%output_station_info
    output_station_mapping=cfg%output_station_mapping
    output_fields_block1=cfg%output_fields_block1
    output_fields_block2=cfg%output_fields_block2
    station_xyz_index=cfg%station_xyz_index
    station_number_in_list=cfg%station_number_in_list
    station_number_in_filename=cfg%station_number_in_filename
    station_use_block_subdirectories=cfg%station_use_block_subdirectories
    station_add_header=cfg%station_add_header
    station_add_metadata=cfg%station_add_metadata
    station_list=cfg%station_list
    station_list_file=cfg%station_list_file
    station_file_directory=cfg%station_file_directory
    station_output_order=cfg%station_output_order
    common_stations_blocks=cfg%common_stations_blocks
    stride_fields=cfg%stride_fields

    S%output_exact_moment = output_exact_moment
    S%output_seismograms = output_seismograms
      S%output_station_info = output_station_info
    S%output_fields_block1 = output_fields_block1
    S%output_fields_block2 = output_fields_block2
    S%stride_fields = stride_fields
      S%station_xyz_index = station_xyz_index
    S%station_name_in_list=cfg%station_name_in_list
    S%station_number_in_list = station_number_in_list
    S%station_number_in_filename = station_number_in_filename
    S%station_use_block_subdirectories = station_use_block_subdirectories
    common_stations_blocks = trim(adjustl(lower_text(common_stations_blocks)))
    if (.not.cfg%split_receivers) then
    select case (trim(common_stations_blocks))
    case ('both')
       S%station_both_blocks = .true.
    case ('block1')
       S%station_both_blocks = .false.
       if (S%block_num /= 1) then
          S%output_seismograms = .false.
          S%output_exact_moment = .false.
       end if
    case ('block2')
       S%station_both_blocks = .false.
       if (S%block_num /= 2) then
          S%output_seismograms = .false.
          S%output_exact_moment = .false.
       end if
    case default
       write(*,'(A,A,A)') 'error: invalid common_stations_blocks "', &
            trim(common_stations_blocks), &
            '"; expected block1, block2, or both'
       error stop 'invalid common_stations_blocks'
    end select
    else
       S%station_both_blocks = cfg%append_block
    endif
    if (station_number_in_filename .and. .not.station_number_in_list) then
       error stop 'station_number_in_filename requires station_number_in_list'
    end if
    call parse_station_output_order(station_output_order, S%station_output_order, stat)
    if (stat /= 0) then
       write(*,'(A,A,A)') 'error: invalid station_output_order "', &
            trim(station_output_order), '"; expected each of t, vx, vy, vz exactly once'
       error stop 'invalid station_output_order'
    end if

    any_station_outputs = S%output_seismograms .or. S%output_exact_moment .or. &
                         S%output_fields_block1 .or. S%output_fields_block2
    block_output_directory = ''
    station_output_directory = ''

    if (any_station_outputs) then
       ! Ensure output directories exist (safe even if called by many MPI ranks)
       if (len_trim(station_file_directory) > 0 .and. trim(adjustl(station_file_directory)) /= '.') then
          call make_station_directory(trim(adjustl(station_file_directory)))
       end if

       write(block_str,'(i0)') S%block_num
       if (len_trim(station_file_directory) > 0 .and. trim(adjustl(station_file_directory)) /= '.') then
          block_output_directory = trim(adjustl(station_file_directory)) // '/block' // trim(adjustl(block_str))
       else
          block_output_directory = 'block' // trim(adjustl(block_str))
       end if
       if (S%output_fields_block1 .or. S%output_fields_block2 .or. &
           ((S%output_seismograms .or. S%output_exact_moment) .and. &
            S%station_use_block_subdirectories)) then
          call make_station_directory(trim(block_output_directory))
       end if

       if (S%station_use_block_subdirectories) then
          station_output_directory = block_output_directory
       else if (len_trim(station_file_directory) > 0) then
          station_output_directory = trim(adjustl(station_file_directory))
       else
          station_output_directory = '.'
       end if
    end if

    ! initialization for seismogram output

    if (S%output_seismograms) then

       ! Determine source of station list (infile or external file)
       if (trim(station_list) == 'extfile' .or. trim(station_list) == 'exfile') then
          ! Read station list from external file
          if (trim(station_list_file) == '') then
             stop 'error: station_list_file must be specified when station_list=extfile'
          end if
          
          open(newunit=external_file_unit, file=trim(station_list_file), status='old', iostat=stat)
          if (stat /= 0) then
             print *, 'error: cannot open station_list_file: ', trim(station_list_file)
             stop
          end if
          
          ! Count number of stations in external file
          stat = 0
          do
             read(external_file_unit,'(a)',iostat=stat) temp
             if (stat /= 0) exit
             if (temp=='!---begin:station_list---') exit
          end do
          
          if (stat == 0 .and. temp=='!---begin:station_list---') then
             S%nstations = 0
             do
                read(external_file_unit,'(a)',end=100) temp
                if (temp=='!---end:station_list---') exit
                S%nstations = S%nstations+1
             end do
          else
             close(external_file_unit)
             stop 'error: !---begin:station_list--- not found in external station file'
          end if
          
       else
          ! Default: Read station list from input file
          ! move to start of station list

          rewind(input)
          stat = 0
          do
             read(input,'(a)',iostat=stat) temp
             if (stat /= 0) exit
             if (temp=='!---begin:station_list---') exit
          end do

          if (stat == 0 .and. temp=='!---begin:station_list---') then
             S%nstations = 0
             do
                read(input,'(a)',end=100) temp
                if (temp=='!---end:station_list---') exit
                S%nstations = S%nstations+1
             end do
          else
             ! Backwards-compatible fallback: per-block lists
             if (S%block_num == 1) then
                rewind(input)
                do
                   read(input,'(a)',end=100) temp
                   if (temp=='!---begin:station_listU---') exit
                end do

                S%nstations = 0
                do
                   read(input,'(a)',end=100) temp
                   if (temp=='!---end:station_listU---') exit
                   S%nstations = S%nstations+1
                end do
             end if

             if (S%block_num == 2) then
                rewind(input)
                do
                   read(input,'(a)',end=100) temp
                   if (temp=='!---begin:station_listV---') exit
                end do

                S%nstations = 0
                do
                   read(input,'(a)',end=100) temp
                   if (temp=='!---end:station_listV---') exit
                   S%nstations = S%nstations+1
                end do
             end if
          end if
       end if

       if (S%output_station_info .and. is_master()) then
          station_info_lines = ''
          station_info_lines(1) = 'Station output information'
          station_info_lines(2) = 'station_list: ' // trim(adjustl(station_list))
          if (trim(adjustl(station_list)) == 'extfile' .or. trim(adjustl(station_list)) == 'exfile') then
             station_info_lines(3) = 'station_list_file: ' // trim(adjustl(station_list_file))
          else
             station_info_lines(3) = 'station_list_file: (from input file)'
          end if
          station_info_lines(4) = 'station_file_directory: ' // trim(adjustl(station_file_directory))
          station_info_lines(5) = 'station_xyz_index: ' // merge('T','F',S%station_xyz_index)
          station_info_lines(6) = 'station_output_order: ' // trim(adjustl(station_output_order))
          station_info_lines(7) = 'station_number_in_list: ' // merge('T','F',S%station_number_in_list)
          station_info_lines(8) = 'station_number_in_filename: ' // &
               merge('T','F',S%station_number_in_filename)
          station_info_lines(9) = 'station_use_block_subdirectories: ' // &
               merge('T','F',S%station_use_block_subdirectories)
          station_info_lines(10) = 'common_stations_blocks: ' // trim(common_stations_blocks)
          if(cfg%split_receivers) station_info_lines(10) = 'interface_stations: '//trim(cfg%interface_stations)
          station_info_lines(11) = 'station_add_header: ' // merge('T','F',station_add_header)
          station_info_lines(12) = 'station_add_metadata: ' // merge('T','F',station_add_metadata)
          station_info_lines(13) = 'output_station_mapping: ' // merge('T','F',output_station_mapping)
          write(station_info_lines(14),'(a,i0)') 'nstations (this block): ', S%nstations
          call boxed_lines(14, station_info_lines(1:14), 78)
       end if

       ! allocate station indices array and output file unit array

       if (S%output_exact_moment) then
         allocate(S%file_unit(2*S%nstations))
       else
         allocate(S%file_unit(S%nstations))
       end if
           
         allocate(S%i(S%nstations),S%j(S%nstations),S%k(S%nstations))
         allocate(S%station_number(S%nstations),S%station_names(S%nstations))
         S%station_names=""
         allocate(S%i_phys(S%nstations),S%j_phys(S%nstations), &
              S%k_phys(S%nstations)) 

       ! read station list again, this time storing station indices

       if (trim(station_list) == 'extfile' .or. trim(station_list) == 'exfile') then
          ! Read from external file
          rewind(external_file_unit)
          stat = 0
          do
             read(external_file_unit,'(a)',iostat=stat) temp
             if (stat /= 0) exit
             if (temp=='!---begin:station_list---') exit
          end do
          
          if (S%nstations > 0) then
             do n = 1,S%nstations
                call read_station_record(external_file_unit, S%station_number_in_list, &
                     n, S%station_number(n), S%i_phys(n), S%j_phys(n), S%k_phys(n), &
                     trim(station_list_file), S%block_num, S%station_name_in_list, S%station_names(n))
             end do
          end if
          
          close(external_file_unit)
          
       else
          ! Read from input file
          rewind(input)

          stat = 0
          do
             read(input,'(a)',iostat=stat) temp
             if (stat /= 0) exit
             if (temp=='!---begin:station_list---') exit
          end do

          if (.not. (stat == 0 .and. temp=='!---begin:station_list---')) then
             if(S%block_num == 1) then
                rewind(input)
                do
                   read(input,'(a)',end=100) temp
                   if (temp=='!---begin:station_listU---') exit
                end do
             end if

             if(S%block_num == 2) then
                rewind(input)
                do
                   read(input,'(a)',end=100) temp
                   if (temp=='!---begin:station_listV---') exit
                end do
             end if
          end if

          if(S%nstations > 0) then
             do n = 1,S%nstations
                call read_station_record(input, S%station_number_in_list, n, &
                     S%station_number(n), S%i_phys(n), S%j_phys(n), S%k_phys(n), &
                     'input station list', S%block_num, S%station_name_in_list, S%station_names(n))
                !print *, S%i_phys(n),S%j_phys(n),S%k_phys(n)
             end do
          end if
       end if

       if(S%station_name_in_list) then
          do n=1,S%nstations
             if(count(lowercase(S%station_names)==lowercase(S%station_names(n)))>1) call fatal_local('CFG-STATION-NAME', &
                  'Duplicate station name: '//trim(S%station_names(n)), 'init_seismogram')
          enddo
       endif
       if (S%station_number_in_filename) then
          do n = 1,S%nstations
             if (count(S%station_number == S%station_number(n)) > 1) then
                write(*,'(A,I0)') 'error: duplicate station number used for output filename: ', &
                     S%station_number(n)
                error stop 'duplicate station number'
             end if
          end do
       end if
       
       
       !xmin = minval(G%X(mx:px,my:py,mz:pz,1))
       !xmax = maxval(G%X(mx:px,my:py,mz:pz,1))
       !ymin = minval(G%X(mx:px,my:py,mz:pz,2))
       !ymax = maxval(G%X(mx:px,my:py,mz:pz,2))
       !zmin = minval(G%X(mx:px,my:py,mz:pz,3))
       !zmax = maxval(G%X(mx:px,my:py,mz:pz,3))
       S%i(:) = -10000
       S%j(:) = -10000
       S%k(:) = -10000
       
       call Find_Coordinates(G%X, S%i_phys, S%j_phys, S%k_phys, &
            S%i, S%j, S%k, S%nstations, mx, my, mz, px, py, pz, &
            output_station_mapping.and..not.cfg%split_receivers)

       if(cfg%split_receivers) then
          do n=1,S%nstations
             ! The two blocks meet at q=nq in block 1 and q=1 in block 2.
             bnd=1
             if(S%block_num==1) bnd=G%C%nq
             local_distance=huge(1.0_wp)
             if(S%i(n)>0) local_distance=sum((G%X(S%i(n),S%j(n),S%k(n),1:3)- &
                  [S%i_phys(n),S%j_phys(n),S%k_phys(n)])**2)
             call MPI_Allreduce(local_distance,best,1,MPI_DOUBLE_PRECISION,MPI_MIN,G%C%comm,ierr)
             candidate=huge(1)
             if(local_distance==best.and.S%i(n)>0) candidate=G%C%rank
             call MPI_Allreduce(candidate,winner,1,MPI_INTEGER,MPI_MIN,G%C%comm,ierr)
             if(G%C%rank/=winner) S%i(n)=-10000
             if(nblocks==2.and.S%i(n)==bnd) then
                if((cfg%interface_stations=='block_1'.and.S%block_num==2).or. &
                   (cfg%interface_stations=='block_2'.and.S%block_num==1)) S%i(n)=-10000
             endif
             if(output_station_mapping.and.S%i(n)>0) then
                write(*,'(A,I0,A,I0,A,3(I0,1X),A,3(ES14.6E3,1X))') &
                     'station ',n,', block ',S%block_num,': selected grid indices=',S%i(n),S%j(n),S%k(n), &
                     ', grid_xyz=',G%X(S%i(n),S%j(n),S%k(n),1:3)
             endif
          enddo
       endif
       allocate(S%station_paths(2*S%nstations))
       S%station_paths=''
       if(cfg%split_receivers) then
          allocate(S%station_requests(S%nstations),S%station_mappings(S%nstations))
          S%station_requests=''; S%station_mappings=''
          do n=1,S%nstations
             write(temp,'(I0,3(1X,ES24.16E3))') S%station_number(n), &
                  S%i_phys(n),S%j_phys(n),S%k_phys(n)
             if(S%station_name_in_list) write(temp,'(A,3(1X,ES24.16E3))') trim(S%station_names(n)), &
                  S%i_phys(n),S%j_phys(n),S%k_phys(n)
             if(G%C%rank==0) S%station_requests(n)=temp
             if(S%i(n)>0) S%station_mappings(n)=temp
          enddo
       endif
       S%station_add_header=station_add_header
       S%station_add_metadata=station_add_metadata
       ! open file units for output

       do n = 1,S%nstations
          if (S%i(n) > 0 .and. S%j(n) > 0 .and. S%k(n) > 0) then
            if(S%station_name_in_list) then
               filename=trim(S%station_names(n))//'.dat'
            else if (S%station_number_in_filename) then
               write(filename,'(a,i0,a)') trim(adjustl(name)) // '_station-', &
                    S%station_number(n), '.dat'
            else if (S%station_xyz_index) then
               write(xs,'(f20.3)') S%i_phys(n)
               write(ys,'(f20.3)') S%j_phys(n)
               write(zs,'(f20.3)') S%k_phys(n)
               write(filename,'(a,a,a,a,a,a)') trim(adjustl(name)) // '_', &
                    trim(adjustl(xs)),'_',trim(adjustl(ys)),'_',trim(adjustl(zs))//'.dat'
            else
               write(filename,'(a,i0,a,i0,a,i0,a,i0,a)') trim(adjustl(name)) // '_', &
                    S%i(n),'_',S%j(n),'_',S%k(n),'_block',S%block_num,'.dat'
            end if
            if(cfg%split_receivers.and..not.(S%station_number_in_filename.or.S%station_xyz_index.or.S%station_name_in_list)) then
               call grid_station_filename(filename,name,S%i(n),S%j(n),S%k(n),cfg%append_block,S%block_num)
            endif
            if (S%station_both_blocks .and. &
                (S%station_number_in_filename .or. S%station_xyz_index .or. S%station_name_in_list)) then
               call append_block_suffix(filename, S%block_num)
            end if
            filepath = trim(station_output_directory) // '/' // trim(filename)
            S%station_paths(n)=filepath
          end if
       end do

       if( S%output_exact_moment) then
       do n = 1,S%nstations
          if (S%i(n) > 0 .and. S%j(n) > 0 .and. S%k(n) > 0) then
            if(S%station_name_in_list) then
               filename=trim(S%station_names(n))//'_exact.dat'
            else if (S%station_number_in_filename) then
               write(filename,'(a,i0,a)') trim(adjustl(name)) // '_exact_station-', &
                    S%station_number(n), '.dat'
            else if (S%station_xyz_index) then
               write(xs,'(f20.3)') S%i_phys(n)
               write(ys,'(f20.3)') S%j_phys(n)
               write(zs,'(f20.3)') S%k_phys(n)
               write(filename,'(a,a,a,a,a,a)') trim(adjustl(name)) // '_exact_', &
                    trim(adjustl(xs)),'_',trim(adjustl(ys)),'_',trim(adjustl(zs))//'.dat'
            else
               write(filename,'(a,i0,a,i0,a,i0,a,i0,a)') trim(adjustl(name)) // '_exact_', &
                    S%i(n),'_',S%j(n),'_',S%k(n),'_block',S%block_num,'.dat'
            end if
            if(cfg%split_receivers.and..not.(S%station_number_in_filename.or.S%station_xyz_index.or.S%station_name_in_list)) then
               call grid_station_filename(filename,trim(name)//"_exact", &
                    S%i(n),S%j(n),S%k(n),cfg%append_block,S%block_num)
            endif
            if (S%station_both_blocks .and. &
                (S%station_number_in_filename .or. S%station_xyz_index .or. S%station_name_in_list)) then
               call append_block_suffix(filename, S%block_num)
            end if
            filepath = trim(station_output_directory) // '/' // trim(filename)
            S%station_paths(S%nstations+n)=filepath
          end if
       end do 
       end if

    end if

    ! initialization for body fields output

    if (S%output_fields_block1.or.S%output_fields_block2) then

       ! open file units for output

       do n = 1,9

          select case(n)
          case(1)
             field = 'vx'
          case(2)
             field = 'vy'
          case(3)
             field = 'vz'
          case(4)
             field = 'sxx'
          case(5)
             field = 'sxy'
          case(6)
             field = 'sxz'
          case(7)
             field = 'syy'
          case(8)
             field = 'syz'
          case(9)
             field = 'szz'
          end select

          if (S%output_fields_block1 .and. S%block_num == 1) then
             write(filename,'(a,a,a)') trim(adjustl(name)) // '_block1_',trim(adjustl(field)),'.dat'
             filepath = trim(block_output_directory) // '/' // trim(filename)
             open(newunit=S%file_unit_block1(n),file=trim(filepath))
          end if

          if (S%output_fields_block2 .and. S%block_num == 2) then
             write(filename,'(a,a,a)') trim(adjustl(name)) // '_block2_',trim(adjustl(field)),'.dat'
             filepath = trim(block_output_directory) // '/' // trim(filename)
             open(newunit=S%file_unit_block2(n),file=trim(filepath))
          end if

       end do

    end if

    return

100 stop 'error reading seismogram station list'

  end subroutine init_seismogram


  subroutine make_station_directory(path)
    use diagnostics, only: fatal_local
    character(*),intent(in)::path
    character(:),allocatable::quoted
    integer::i,exit_status,command_status
    character(512)::message
    ! Single-quote the path so shell substitutions in input cannot execute.
    quoted=achar(39)
    do i=1,len_trim(path)
      if(path(i:i)==achar(39)) then
        quoted=quoted//achar(39)//achar(34)//achar(39)//achar(34)//achar(39)
      else
        quoted=quoted//path(i:i)
      endif
    enddo
    quoted=quoted//achar(39)
    call execute_command_line('mkdir -p -- '//quoted,exitstat=exit_status, &
         cmdstat=command_status,cmdmsg=message)
    if(command_status/=0) call fatal_local('CFG-RECEIVER-DIRECTORY',trim(message),'make_station_directory')
    if(exit_status/=0) call fatal_local('CFG-RECEIVER-DIRECTORY', &
         'Cannot create output directory: '//path,'make_station_directory')
  end subroutine

  subroutine grid_station_filename(filename,name,i,j,k,append_block,block)
    character(*),intent(out)::filename
    character(*),intent(in)::name
    integer,intent(in)::i,j,k,block
    logical,intent(in)::append_block
    write(filename,'(A,I0,A,I0,A,I0,A)') trim(name)//'_',i,'_',j,'_',k,'.dat'
    if(append_block) call append_block_suffix(filename,block)
  end subroutine

  subroutine open_station_files(S,G)
    use diagnostics, only: fatal_local
    type(seismogram_type),intent(inout)::S
    type(block_grid_t),intent(in)::G
    integer::n,m,stat
    character(512)::message
    if(.not.allocated(S%station_paths)) return
    do m=1,size(S%station_paths)
      if(len_trim(S%station_paths(m))==0) cycle
      n=mod(m-1,S%nstations)+1
      open(newunit=S%file_unit(m),file=trim(S%station_paths(m)),status='replace',iostat=stat,iomsg=message)
      if(stat/=0) call fatal_local('CFG-RECEIVER-OPEN',trim(message),'open_station_files')
      call write_station_preamble(S%file_unit(m),S%station_add_header,S%station_add_metadata, &
           S%station_output_order,S%station_number_in_list,S%station_number(n), &
           S%i_phys(n),S%j_phys(n),S%k_phys(n),S%i(n),S%j(n),S%k(n), &
           G%X(S%i(n),S%j(n),S%k(n),1:3),S%station_names(n))
    enddo
  end subroutine

  subroutine check_station_paths(stations)
    use diagnostics, only: fatal_local
    type(seismogram_type),intent(in)::stations(:)
    character(513),allocatable::local(:),allpaths(:)
    integer,allocatable::counts(:),offsets(:)
    integer::i,j,k,n,total,ranks,ierr
    logical::found
    n=0
    do i=1,size(stations)
      if(allocated(stations(i)%station_paths)) n=n+count(stations(i)%station_paths/='')
      if(allocated(stations(i)%station_requests)) n=n+count(stations(i)%station_requests/='')+ &
          count(stations(i)%station_mappings/='')
    enddo
    allocate(local(n)); k=0
    do i=1,size(stations)
      if(.not.allocated(stations(i)%station_paths)) cycle
      do j=1,size(stations(i)%station_paths)
        if(stations(i)%station_paths(j)=='') cycle
        k=k+1; local(k)='P'//trim(stations(i)%station_paths(j))
      enddo
      if(allocated(stations(i)%station_requests)) then
        do j=1,size(stations(i)%station_requests)
          if(stations(i)%station_requests(j)/='') then
            k=k+1; local(k)='R'//trim(stations(i)%station_requests(j))
          endif
          if(stations(i)%station_mappings(j)/='') then
            k=k+1; local(k)='M'//trim(stations(i)%station_mappings(j))
          endif
        enddo
      endif
    enddo
    call MPI_Comm_size(MPI_COMM_WORLD,ranks,ierr)
    allocate(counts(ranks),offsets(ranks))
    n=n*513
    call MPI_Allgather(n,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    offsets(1)=0
    do i=2,ranks
      offsets(i)=offsets(i-1)+counts(i-1)
    enddo
    total=sum(counts)/513
    allocate(allpaths(total))
    call MPI_Allgatherv(local,n,MPI_CHARACTER,allpaths,counts,offsets,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
    do i=1,total
      if(allpaths(i)(1:1)=='R') then
        found=.false.
        do j=1,total
          if(allpaths(j)(1:1)=='M'.and.allpaths(i)(2:)==allpaths(j)(2:)) found=.true.
        enddo
        if(.not.found) call fatal_local('CFG-RECEIVER-UNMAPPED', &
          'No selected mesh node for station (number x y z): '//trim(allpaths(i)(2:)), 'check_station_paths')
      endif
      if(allpaths(i)(1:1)/='P') cycle
      do j=i+1,total
        if(lowercase(allpaths(i))==lowercase(allpaths(j))) call fatal_local('CFG-RECEIVER-COLLISION', &
          'Duplicate station output path: '//trim(allpaths(i)(2:)),'check_station_paths')
      enddo
    enddo
  end subroutine

  subroutine write_seismogram(S, t, F)

    use mpi3dbasic, only : rank
    use diagnostics, only : fatal_local
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none

    type(seismogram_type),intent(in) :: S
    real(kind = wp),intent(in) :: t
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(in) :: F
    integer :: n, bad_component
    real(kind = wp) :: station_values(4)
    character(len=256) :: message
    character(len=4), parameter :: component_name(4) = [character(len=4) :: &
         'time', 'vx', 'vy', 'vz']

    if (S%output_seismograms .and. S%nstations > 0) then
 
       do n = 1,S%nstations
          if (S%i(n) > 0 .and. S%j(n) > 0 .and. S%k(n) > 0) then
            station_values = [t, F(S%i(n),S%j(n),S%k(n),1:3)]
            if (.not.all(ieee_is_finite(station_values))) then
               bad_component = findloc(.not.ieee_is_finite(station_values), .true., dim=1)
               write(message,'(A,I0,A,I0,A,I0,A,I0,A,A,A,ES24.16E3)') &
                    'Non-finite station output at station ', n, ' (i,j,k)=(', &
                    S%i(n), ',', S%j(n), ',', S%k(n), '), component ', &
                    trim(component_name(bad_component)), ', value=', station_values(bad_component)
               call fatal_local('RUN-STATION-001', trim(message), &
                    'write_seismogram', S%block_num)
            end if
            write(S%file_unit(n),'(ES24.16E3,3(1X,ES24.16E3))') &
                 station_values(S%station_output_order)
          end if
       end do

    end if

    if (S%output_fields_block1 .and. S%block_num == 1) then
       do n = 1,9
          write(S%file_unit_block1(n),*) F(:,:,:,n)
          print *, S%file_unit_block1(n)
          flush(S%file_unit_block1(n))
       end do

    end if

  end subroutine write_seismogram

  subroutine parse_station_output_order(order_text, order, stat)

    implicit none

    character(*), intent(in) :: order_text
    integer, intent(out) :: order(4)
    integer, intent(out) :: stat
    integer :: i, n, component
    logical :: seen(4)
    character :: first, second

    order = 0
    seen = .false.
    stat = 0
    i = 1
    n = 0

    do while (i <= len_trim(order_text))
       first = lower_ascii(order_text(i:i))
       if (first == ' ' .or. first == ',' .or. first == char(9)) then
          i = i + 1
          cycle
       end if

       if (first == 't') then
          component = 1
          i = i + 1
       else if (first == 'v' .and. i < len_trim(order_text)) then
          second = lower_ascii(order_text(i+1:i+1))
          select case (second)
          case ('x')
             component = 2
          case ('y')
             component = 3
          case ('z')
             component = 4
          case default
             stat = 1
             return
          end select
          i = i + 2
       else
          stat = 1
          return
       end if

       n = n + 1
       if (n > 4 .or. seen(component)) then
          stat = 1
          return
       end if
       order(n) = component
       seen(component) = .true.
    end do

    if (n /= 4 .or. .not.all(seen)) stat = 1

  end subroutine parse_station_output_order

  subroutine append_block_suffix(filename, block_num)

    implicit none

    character(*), intent(inout) :: filename
    integer, intent(in) :: block_num
    integer :: extension_position
    character(32) :: suffix

    extension_position = index(trim(filename), '.dat', back=.true.)
    if (extension_position == 0) then
       error stop 'station filename is missing .dat extension'
    end if
    write(suffix,'(A,I0)') '_block', block_num
    filename = filename(:extension_position-1) // trim(suffix) // '.dat'

  end subroutine append_block_suffix

  subroutine write_station_preamble(file_unit, add_header, add_metadata, order, &
       has_station_number, station_number, x, y, z, i, j, k, grid_xyz, station_name)

    implicit none

    integer, intent(in) :: file_unit, order(4), station_number, i, j, k
    logical, intent(in) :: add_header, add_metadata, has_station_number
    real(kind=wp), intent(in) :: x, y, z
    character(*),intent(in) :: station_name
    real(kind=wp), intent(in) :: grid_xyz(3)
    real(kind=wp) :: mapping_distance
    integer :: n
    character(32), parameter :: component_name(4) = [character(32) :: &
         't', 'vx', 'vy', 'vz']
    character(256) :: header

    if (add_metadata) then
       if(len_trim(station_name)>0) write(file_unit,'(A)') "# station_name: "//trim(station_name)
       if (has_station_number) write(file_unit,'(A,I0)') '# station_number: ', station_number
       write(file_unit,'(A,3(1X,ES24.16E3))') '# x y z:', x, y, z
       mapping_distance = sqrt((grid_xyz(1)-x)**2 + (grid_xyz(2)-y)**2 + &
            (grid_xyz(3)-z)**2)
       write(file_unit,'(A,3(1X,I0))') '# grid_i j k:', i, j, k
       write(file_unit,'(A,3(1X,ES24.16E3))') '# grid_x y z:', grid_xyz
       write(file_unit,'(A,1X,ES24.16E3)') '# mapping_distance:', mapping_distance
    end if

    if (add_header) then
       header = '#'
       do n = 1,4
          header = trim(header) // ' ' // trim(component_name(order(n)))
       end do
       write(file_unit,'(A)') trim(header)
    end if

  end subroutine write_station_preamble

  subroutine read_station_record(file_unit, has_station_number, row_number, &
       station_number, x, y, z, source, block_num, has_station_name, station_name)

    use diagnostics, only : fatal_local
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    implicit none

    integer, intent(in) :: file_unit, row_number, block_num
    logical, intent(in) :: has_station_number, has_station_name
    character(*),intent(out)::station_name
    character(256)::name_token
    integer, intent(out) :: station_number
    real(kind=wp), intent(out) :: x, y, z
    character(*), intent(in) :: source
    integer :: stat
    character(256) :: line
    character(512) :: iomsg, message

    station_name=''; name_token=''
    if(has_station_name) then
       x=ieee_value(0.0_wp,ieee_quiet_nan); y=x; z=x
    endif
    read(file_unit,'(A)',iostat=stat,iomsg=iomsg) line
    if (stat == 0) then
       if(has_station_name) then
          station_number=row_number
          read(line,*,iostat=stat,iomsg=iomsg) name_token,x,y,z
          if(stat==0) then
             if(.not.valid_station_name(name_token)) call fatal_local('CFG-STATION-NAME', &
                  'Station names must be 1-64 letters, digits, hyphens or underscores', 'read_station_record')
             station_name=name_token
             if(.not.all(ieee_is_finite([x,y,z]))) stat=1
          endif
       else if (has_station_number) then
          read(line,*,iostat=stat,iomsg=iomsg) station_number, x, y, z
       else
          station_number = row_number
          read(line,*,iostat=stat,iomsg=iomsg) x, y, z
       end if
    end if

    if (stat /= 0) then
       if(has_station_name) then
          write(message,'(A,I0,A,A)') 'Cannot parse named station row ',row_number,'; expected name x y z: ',trim(line)
       else if (has_station_number) then
          write(message,'(A,I0,A,A,A,A)') 'Cannot parse station row ', row_number, &
               ' from ', trim(source), '; expected: integer_station_number x y z. Row: ', &
               trim(line)
       else
          write(message,'(A,I0,A,A,A,A)') 'Cannot parse station row ', row_number, &
               ' from ', trim(source), '; expected: x y z. Row: ', trim(line)
       end if
       call fatal_local('CFG-STATION-001', trim(message), 'read_station_record', &
            block_num, root_only=.true.)
    end if

  end subroutine read_station_record

  pure function lower_text(value) result(lower)

    implicit none

    character(*), intent(in) :: value
    character(len(value)) :: lower
    integer :: i

    do i = 1,len(value)
       lower(i:i) = lower_ascii(value(i:i))
    end do

  end function lower_text

  pure function lower_ascii(value) result(lower)

    implicit none

    character, intent(in) :: value
    character :: lower
    integer :: code

    code = iachar(value)
    if (code >= iachar('A') .and. code <= iachar('Z')) then
       lower = achar(code + iachar('a') - iachar('A'))
    else
       lower = value
    end if

  end function lower_ascii

  subroutine destroy_seismogram(S)

    implicit none

    type(seismogram_type),intent(inout) :: S

    integer :: n

    if (S%output_seismograms) then

      if(allocated(S%station_paths)) then
         do n=1,size(S%station_paths)
            if(S%station_paths(n)/='') close(S%file_unit(n))
         enddo
      endif

      if(allocated(S%station_names)) deallocate(S%station_names)
      if(allocated(S%station_paths)) deallocate(S%station_paths)
      if(allocated(S%station_requests)) deallocate(S%station_requests,S%station_mappings)
      if (allocated(S%i)) deallocate(S%i,S%j,S%k,S%file_unit, &
           S%i_phys,S%j_phys,S%k_phys,S%station_number)

    end if

    if (S%output_fields_block1) then

       do n = 1,9
          close(S%file_unit_block1(n))
       end do

    end if

    if (S%output_fields_block2) then

       do n = 1,9
          close(S%file_unit_block2(n))
       end do

    end if

  end subroutine destroy_seismogram



  subroutine Find_Coordinates(XX, x1, y1, z1, x_i, y_j, z_k, nstations, &
       mx, my, mz, px, py, pz, output_station_mapping)

    ! Given the physical positions of the receivers x1, y1, z1
    ! Find the corresponding indices in x_i, y_j, z_k in the mesh XX

    implicit none

    integer, intent(in) :: mx, my, mz, px, py, pz               ! size of the grid-block
    integer, intent(in) :: nstations
    logical, intent(in), optional :: output_station_mapping
    real(kind = wp), dimension(:), intent(in) :: x1, y1, z1                ! receiver  positions
    integer, dimension(:), intent(out) :: x_i, y_j, z_k         ! spatial indices of receiver  positions to be found
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(in) :: XX                  ! grid
    integer :: i, j, k, c,  i0, j0, k0
    real(kind = wp) :: vec(3), dist, mindist,xmin, xmax, ymin, ymax, zmin, zmax
    real(kind = wp) :: hx,hy,hz, dist0
    logical :: print_station_mapping

    print_station_mapping = .true.
    if (present(output_station_mapping)) print_station_mapping = output_station_mapping

    xmin = minval(XX(mx:px,my:py,mz:pz,1))
    xmax = maxval(XX(mx:px,my:py,mz:pz,1))
    ymin = minval(XX(mx:px,my:py,mz:pz,2))
    ymax = maxval(XX(mx:px,my:py,mz:pz,2))
    zmin = minval(XX(mx:px,my:py,mz:pz,3))
    zmax = maxval(XX(mx:px,my:py,mz:pz,3))
    hx = 1d0
    hy = 1d0
    hz = 1d0

    !print *,  xmin, xmax, ymin, ymax, zmin, zmax
    !print *, nstations
    do c = 1, nstations
       
        

       i0 = -9999
       j0 = -9999
       k0 = -9999
          
       mindist = 1.0e8_wp


       !if ((xmin <= x1(c) .and. x1(c) <= xmax) .and. &
       !     (ymin <= y1(c) .and. y1(c) <= ymax) .and. &
       !     (zmin <= z1(c) .and. z1(c) <= zmax)) then

          


          k_loop: do k = mz, pz
             do j = my, py
                do i = mx, px

                   hx = XX(px, j, k, 1)-XX(px-1, j, k, 1)
                   hy = XX(i, py, k, 2)-XX(i, py-1, k, 2)
                   hz = XX(i, j, pz, 3)-XX(i, j, pz-1, 3)

                   if (i<px) hx = XX(i+1, j, k, 1)-XX(i, j, k, 1)
                   if (j<py) hy = XX(i, j+1, k, 2)-XX(i, j, k, 2)
                   if (k<pz) hz = XX(i, j, k+1, 3)-XX(i, j, k, 3)


                   vec = [XX(i, j, k, 1)-x1(c), XX(i, j, k, 2)-y1(c), XX(i, j, k, 3)-z1(c)]
                   dist = sqrt(dot_product(vec, vec))
                   dist0 = sqrt(hx**2 + hy**2 + hz**2)
                   
                   !if (((abs(XX(i, j, k, 1)-x1(c))<=0.8_wp*hx) .and. &
                   !     (abs(XX(i, j, k, 2)-y1(c))<=0.8_wp*hy) .and. &
                   !     (abs(XX(i, j, k, 3)-z1(c))<=0.8_wp*hz)) .and. &
                   !      (dist <= mindist)) then

                   if ((dist <= mindist) .and. (dist <= 0.5*dist0))  then
                   !if ((dist <= mindist))  then

                      i0 = i
                      j0 = j
                      k0 = k

                      x_i(c) = i0
                      y_j(c) = j0
                      z_k(c) = k0

                      mindist = dist

                   end if
                   
                   
                end do
             end do

          end do k_loop
          !print*, mindist, c, i0, j0, k0, XX(i0, j0, k0, 1), XX(i0, j0, k0, 2), XX(i0, j0, k0, 3), x1(c), y1(c), z1(c)
       !end if
       if (print_station_mapping .and. i0 > 0 .and. j0>0 .and. k0 > 0) then
          write(*,'(A,I0,A,ES14.6E3,A,3(I0,1X),A,3(ES14.6E3,1X),A,3(ES14.6E3,1X),A)') &
               'station ', c, ': distance=', mindist, ', indices=(', i0, j0, k0, &
               '), grid_xyz=(', XX(i0,j0,k0,1), XX(i0,j0,k0,2), XX(i0,j0,k0,3), &
               '), requested_xyz=(', x1(c), y1(c), z1(c), ');'
       end if
       
        !x_i(c) = i0
        !y_j(c) = j0
        !z_k(c) = k0
     end do


  end subroutine Find_Coordinates

    subroutine Find_Coordinates_moment(XX, x1, y1, z1, x_i, y_j, z_k, nstations,nq,nr,ns, mx, my, mz, px, py, pz)

    ! Given the physical positions of the receivers x1, y1, z1
    ! Find the corresponding indices in x_i, y_j, z_k in the mesh XX

    implicit none

    integer, intent(in) :: mx, my, mz, px, py, pz, nq, nr, ns               ! size of the grid-block
    integer, intent(in) :: nstations
    real(kind = wp), dimension(:), intent(in) :: x1, y1, z1                ! receiver  positions
    integer, dimension(:), intent(out) :: x_i, y_j, z_k         ! spatial indices of receiver  positions to be found
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(in) :: XX                  ! grid
    integer :: i, j, k, c
    real(kind = wp) :: vec(3), dist, mindist,xmin, xmax, ymin, ymax, zmin, zmax, tolerance, scale
    real(kind = wp) :: hx,hy,hz

    xmin = minval(XX(mx:px,my:py,mz:pz,1))
    xmax = maxval(XX(mx:px,my:py,mz:pz,1))
    ymin = minval(XX(mx:px,my:py,mz:pz,2))
    ymax = maxval(XX(mx:px,my:py,mz:pz,2))
    zmin = minval(XX(mx:px,my:py,mz:pz,3))
    zmax = maxval(XX(mx:px,my:py,mz:pz,3))
    scale = max(1.0_wp, abs(xmin), abs(xmax), abs(ymin), abs(ymax), abs(zmin), abs(zmax))
    tolerance = 1000.0_wp*epsilon(1.0_wp)*scale

    hx = 1d0
    hy = 1d0
    hz = 1d0

    do c = 1, nstations
       x_i(c) = -9999
       y_j(c) = -9999
       z_k(c) = -9999



       if ((xmin-tolerance <= x1(c) .and. x1(c) <= xmax+tolerance) .and. &
            (ymin-tolerance <= y1(c) .and. y1(c) <= ymax+tolerance) .and. &
            (zmin-tolerance <= z1(c) .and. z1(c) <= zmax+tolerance)) then
          
          mindist = 1.0e8_wp
          !        hx = XX(px, j, k, 1)-XX(px-1, j, k, 1)
          !        hy = XX(i, py, k, 1)-XX(i, py-1, k, 2)
          !        hz = XX(i, j, pz, 1)-XX(i, j, pz-1, 3)


          k_loop: do k = mz, pz
             do j = my, py
                do i = mx, px
                   
                   !print*, i,j,k,XX(i,j,k,1)
                   
                   if (i<nq) hx = XX(i+1, j, k, 1)-XX(i, j, k, 1)
                   if (j<nr) hy = XX(i, j+1, k, 2)-XX(i, j, k, 2)
                   if (k<ns) hz = XX(i, j, k+1, 3)-XX(i, j, k, 3)
                   
                   vec = [XX(i, j, k, 1)-x1(c), XX(i, j, k, 2)-y1(c), XX(i, j, k, 3)-z1(c)]
                   dist = sqrt(dot_product(vec, vec))
                   
                   !                    if (((abs(XX(i, j, k, 1)-x1(c))<=0.9d0*hx) .and. &
                   !                         (abs(XX(i, j, k, 2)-y1(c))<=0.9d0*hy) .and. &
                   !                         (abs(XX(i, j, k, 3)-z1(c))<=0.9d0*hz))) then
                   
                   if (mx == 1 .and. px == 25 .and. my == 1 .and. py == 25) then
                      !print*,XX(i,j,k,1)
                      ! print*,'dist_min = :',i,j,k,dist
                   end if
                   if (dist <= mindist) then
                      
                      x_i(c) = i
                      y_j(c) = j
                      z_k(c) = k
                      
                      mindist = dist
                      
                   end if
                end do
             end do
          end do k_loop
       end if
    end do
       
  end subroutine Find_Coordinates_moment

end module seismogram
