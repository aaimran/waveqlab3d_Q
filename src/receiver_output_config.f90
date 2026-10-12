module receiver_output_config
  use mpi
  implicit none
  type :: receiver_config_t
    logical :: output_exact_moment = .false.
    logical :: output_seismograms = .false.
    logical :: output_station_info = .true.
    logical :: output_station_mapping = .true.
    logical :: output_fault_topo = .false.
    logical :: output_fields_block1 = .false.
    logical :: output_fields_block2 = .false.
    logical :: station_xyz_index = .false.
    logical :: station_name_in_list = .false.
    logical :: station_number_in_list = .false.
    logical :: station_number_in_filename = .false.
    logical :: station_use_block_subdirectories = .true.
    logical :: station_add_header = .false.
    logical :: station_add_metadata = .false.
    logical :: append_block = .true.
    logical :: split_receivers = .false.
    character(256) :: station_list = 'infile'
    character(256) :: station_list_file = ''
    character(256) :: station_file_directory = 'seismogram'
    character(256) :: station_output_order = 't vx vy vz'
    character(256) :: common_stations_blocks = 'both'
    character(256) :: interface_stations = 'both'
    integer :: stride_fields = 1
  end type
contains
  pure elemental function lowercase(text) result(out)
    character(*), intent(in) :: text
    character(len(text)) :: out
    integer :: i, c
    out=text
    do i=1,len(text)
      c=iachar(text(i:i))
      if(c>=65.and.c<=90) out(i:i)=achar(c+32)
    end do
  end function

  subroutine read_receiver_config(unit, cfg, stat, message)
    integer,intent(in) :: unit
    type(receiver_config_t),intent(out) :: cfg
    integer,intent(out) :: stat
    character(*),intent(out) :: message
    character(4096) :: line, clean
    character(:),allocatable :: legacy_text
    integer :: ios, section, old_count, new_count
    logical :: &
       output_exact_moment, output_seismograms, output_station_info, &
       output_station_mapping, output_fault_topo, output_fields_block1, &
       output_fields_block2, station_xyz_index, station_number_in_list, &
       station_number_in_filename, station_use_block_subdirectories, station_add_header, &
       station_add_metadata, append_block, split_receivers, station_name_in_list
    character(256) :: &
       station_list, station_list_file, station_file_directory, &
       station_output_order, common_stations_blocks, interface_stations
    integer :: stride_fields
    namelist /output_list/ &
       output_exact_moment, output_seismograms, output_station_info, &
       output_station_mapping, output_fault_topo, output_fields_block1, &
       output_fields_block2, station_xyz_index, station_number_in_list, &
       station_number_in_filename, station_use_block_subdirectories, station_add_header, &
       station_add_metadata, station_list, station_list_file, &
       station_file_directory, station_output_order, common_stations_blocks, &
       stride_fields
    namelist /output_receiver_stations/ &
       output_seismograms, output_station_info, output_station_mapping, &
       station_xyz_index, station_number_in_list, station_number_in_filename, &
       station_use_block_subdirectories, station_add_header, station_add_metadata, &
       station_list, station_list_file, station_file_directory, &
       station_output_order, append_block, interface_stations, station_name_in_list
    cfg=receiver_config_t()
    stat=0; message=""
    output_exact_moment=cfg%output_exact_moment
    output_seismograms=cfg%output_seismograms
    output_station_info=cfg%output_station_info
    output_station_mapping=cfg%output_station_mapping
    output_fault_topo=cfg%output_fault_topo
    output_fields_block1=cfg%output_fields_block1
    output_fields_block2=cfg%output_fields_block2
    station_xyz_index=cfg%station_xyz_index
    station_name_in_list=cfg%station_name_in_list
    station_number_in_list=cfg%station_number_in_list
    station_number_in_filename=cfg%station_number_in_filename
    station_use_block_subdirectories=cfg%station_use_block_subdirectories
    station_add_header=cfg%station_add_header
    station_add_metadata=cfg%station_add_metadata
    append_block=cfg%append_block
    split_receivers=cfg%split_receivers
    station_list=cfg%station_list
    station_list_file=cfg%station_list_file
    station_file_directory=cfg%station_file_directory
    station_output_order=cfg%station_output_order
    common_stations_blocks=cfg%common_stations_blocks
    interface_stations=cfg%interface_stations
    stride_fields=cfg%stride_fields
    section=0; old_count=0; new_count=0; legacy_text=''
    rewind(unit)
    do
      read(unit,'(A)',iostat=ios) line
      if(ios/=0) exit
      clean=trim(adjustl(lowercase(line)))
      ! Strip comments only outside quoted strings.
      call remove_comment(clean)
      if(index(clean,'&output_receiver_stations')==1) then
        new_count=new_count+1; section=2
      else if(index(clean,'&output_list')==1) then
        old_count=old_count+1; section=1
      endif
      line=clean
      call remove_quoted_values(line)
      if(section==1) legacy_text=legacy_text//' '//trim(line)
      call remove_quoted_values(clean)
      if(index(clean,'/')>0.or.index(clean,'&end')>0) section=0
    enddo
    if(old_count/=1.or.new_count>1) then
      stat=1; message='Expected one &output_list and at most one &output_receiver_stations'; return
    endif
    rewind(unit)
    read(unit,nml=output_list,iostat=stat,iomsg=message)
    if(stat/=0) return
    if(new_count==1) then
      if(has_assignment(legacy_text,'output_seismograms')) then
        stat=1; message='Receiver setting output_seismograms occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'output_station_info')) then
        stat=1; message='Receiver setting output_station_info occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'output_station_mapping')) then
        stat=1; message='Receiver setting output_station_mapping occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_xyz_index')) then
        stat=1; message='Receiver setting station_xyz_index occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_number_in_list')) then
        stat=1; message='Receiver setting station_number_in_list occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_number_in_filename')) then
        stat=1; message='Receiver setting station_number_in_filename occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_use_block_subdirectories')) then
        stat=1
        message='Receiver setting station_use_block_subdirectories occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_add_header')) then
        stat=1; message='Receiver setting station_add_header occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_add_metadata')) then
        stat=1; message='Receiver setting station_add_metadata occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_list')) then
        stat=1; message='Receiver setting station_list occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_list_file')) then
        stat=1; message='Receiver setting station_list_file occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_file_directory')) then
        stat=1; message='Receiver setting station_file_directory occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'station_output_order')) then
        stat=1; message='Receiver setting station_output_order occurs in &output_list with split receiver section'; return
      endif
      if(has_assignment(legacy_text,'common_stations_blocks')) then
        stat=1; message='Receiver setting common_stations_blocks occurs in &output_list with split receiver section'; return
      endif
      rewind(unit)
      read(unit,nml=output_receiver_stations,iostat=stat,iomsg=message)
      if(stat/=0) return
      split_receivers=.true.
    endif
    station_list=trim(adjustl(lowercase(station_list)))
    interface_stations=trim(adjustl(lowercase(interface_stations)))
    common_stations_blocks=trim(adjustl(lowercase(common_stations_blocks)))
    select case(trim(interface_stations))
    case('block_1','block_2','both')
    case default
      stat=1; message='interface_stations must be block_1, block_2, or both'; return
    end select
    select case(trim(common_stations_blocks))
    case('block1','block2','both')
    case default
      stat=1; message='common_stations_blocks must be block1, block2, or both'; return
    end select
    if(station_name_in_list.and.(station_number_in_list.or.station_number_in_filename.or.station_xyz_index)) then
      stat=1; message='station_name_in_list requires numbered and coordinate filename flags to be F'; return
    endif
    if(station_number_in_filename.and..not.station_number_in_list) then
      stat=1; message='station_number_in_filename requires station_number_in_list'; return
    endif
    cfg%output_exact_moment=output_exact_moment
    cfg%output_seismograms=output_seismograms
    cfg%output_station_info=output_station_info
    cfg%output_station_mapping=output_station_mapping
    cfg%output_fault_topo=output_fault_topo
    cfg%output_fields_block1=output_fields_block1
    cfg%output_fields_block2=output_fields_block2
    cfg%station_xyz_index=station_xyz_index
    cfg%station_name_in_list=station_name_in_list
    cfg%station_number_in_list=station_number_in_list
    cfg%station_number_in_filename=station_number_in_filename
    cfg%station_use_block_subdirectories=station_use_block_subdirectories
    cfg%station_add_header=station_add_header
    cfg%station_add_metadata=station_add_metadata
    cfg%append_block=append_block
    cfg%split_receivers=split_receivers
    cfg%station_list=station_list
    cfg%station_list_file=station_list_file
    cfg%station_file_directory=station_file_directory
    cfg%station_output_order=station_output_order
    cfg%common_stations_blocks=common_stations_blocks
    cfg%interface_stations=interface_stations
    cfg%stride_fields=stride_fields
    if(station_list/='infile'.and.station_list/='extfile'.and.station_list/='exfile') then
      stat=1; message='station_list must be infile, extfile, or exfile'; return
    endif
    if(stride_fields<1) then
      stat=1; message='stride_fields must be positive'; return
    endif
    message=''
  end subroutine

  subroutine remove_comment(line)
    character(*),intent(inout)::line
    character :: quote, c
    integer :: i
    quote=' '
    do i=1,len_trim(line)
      c=line(i:i)
      if(quote/=' ') then
        if(c==quote) quote=' '
      else if(c=="'".or.c=='"') then
        quote=c
      else if(c=='!') then
        line(i:)=' '; exit
      endif
    enddo
  end subroutine

  subroutine remove_quoted_values(line)
    character(*),intent(inout)::line
    character::quote,c
    integer::i
    quote=' '
    do i=1,len_trim(line)
      c=line(i:i)
      if(quote/=' ') then
        line(i:i)=' '
        if(c==quote) quote=' '
      else if(c=="'".or.c=='"') then
        quote=c; line(i:i)=' '
      endif
    enddo
  end subroutine

  logical function has_assignment(text,key) result(found)
    character(*),intent(in)::text,key
    integer :: i,j,c
    found=.false.
    do i=1,len(text)-len(key)+1
      if(text(i:i+len(key)-1)/=key) cycle
      if(i>1) then
        c=iachar(text(i-1:i-1))
        if((c>=97.and.c<=122).or.(c>=48.and.c<=57).or.c==95) cycle
      endif
      j=i+len(key)
      do while(j<=len(text))
        if(text(j:j)/=' '.and.text(j:j)/=achar(9)) exit
        j=j+1
      enddo
      if(j<=len(text)) then
        if(text(j:j)=='=') found=.true.
      endif
    enddo
  end function

  pure logical function valid_station_name(name) result(valid)
    character(*),intent(in)::name
    integer::i,c
    valid=.false.
    if(len_trim(name)<1.or.len_trim(name)>64) return
    do i=1,len_trim(name)
      c=iachar(name(i:i))
      if((c>=65.and.c<=90).or.(c>=97.and.c<=122).or.(c>=48.and.c<=57).or.c==45.or.c==95) cycle
      return
    enddo
    valid=.true.
  end function

  subroutine broadcast_receiver_config(cfg)
    type(receiver_config_t),intent(inout)::cfg
    integer::ierr
    call MPI_Bcast(cfg%output_exact_moment,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%output_seismograms,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%output_station_info,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%output_station_mapping,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%output_fault_topo,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%output_fields_block1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%output_fields_block2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_xyz_index,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_name_in_list,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_number_in_list,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_number_in_filename,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_use_block_subdirectories,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_add_header,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_add_metadata,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%append_block,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%split_receivers,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_list,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_list_file,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_file_directory,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%station_output_order,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%common_stations_blocks,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%interface_stations,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(cfg%stride_fields,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  end subroutine
end module
