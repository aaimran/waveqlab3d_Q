module anelastic_cg_t_comm
 use common,only:wp
 use anelastic_cg_t_types
 use anelastic_cg8_layout,only:cg8_cell_eligible,cg8_mechanism
 use mpi3dcomm,only:cartesian3d_t
 use diagnostics,only:fatal_local
 use mpi
 implicit none
 private
 public::init_cgt_layout,exchange_cgt_stage,cgt_fail
contains
 subroutine cgt_fail(message)
  character(*),intent(in)::message
  call fatal_local('RUN-CGT',message,'projected-cell attenuation')
 end subroutine
 subroutine mpi_check(ierr)
  integer,intent(in)::ierr
  if(ierr/=MPI_SUCCESS) call cgt_fail('MPI cell communication failed')
 end subroutine
 integer function node_rank(point,ranges,me) result(r)
  integer,intent(in)::point(3),ranges(:,0:),me
  integer::i
  r=me
  if(all(point>=ranges(1:3,me)).and.all(point<=ranges(4:6,me))) return
  do i=0,size(ranges,2)-1
    if(all(point>=ranges(1:3,i)).and.all(point<=ranges(4:6,i))) then
      r=i;return
    endif
  enddo
  call cgt_fail('Cell node has no owning rank')
 end function
 pure function slot_point(first,slot) result(point)
  integer,intent(in)::first(3),slot
  integer::point(3),i
  i=slot-1;point=first+[modulo(i,2),modulo(i/2,2),i/4]
 end function
 subroutine init_cgt_layout(s,c)
  type(cgt_state),intent(inout)::s
  type(cartesian3d_t),intent(in)::c
  type(cgt_peer),allocatable::all_peers(:)
  integer,allocatable::ranges(:,:),map(:,:,:),counts(:,:),positions(:,:)
  logical,allocatable::touched(:)
  integer::np,ierr,lo(3),hi(3),anchor(3),last(3),extent(3),first(3),point(3),ijk(3)
  integer::x,y,z,nlocal,nowned,nfull,j,k,r,p,npeer,a,b,d,e,remote_slots,global_remote,max_peers
  integer(kind=8)::local_counts(3),global_counts(3)
  call MPI_Comm_size(c%comm,np,ierr);call mpi_check(ierr)
  s%comm=c%comm;s%rank=c%rank;lo=[c%mq,c%mr,c%ms];hi=[c%pq,c%pr,c%ps]
  allocate(ranges(6,0:np-1),touched(0:np-1),counts(4,0:np-1),positions(4,0:np-1))
  call MPI_Allgather([lo,hi],6,MPI_INTEGER,ranges,6,MPI_INTEGER,c%comm,ierr);call mpi_check(ierr)
  anchor=lo-modulo(modulo(lo,2)-modulo(s%origin,2),2)
  last=hi-modulo(modulo(hi,2)-modulo(s%origin,2),2);extent=(last-anchor)/2+1
  allocate(map(extent(1),extent(2),extent(3)));map=0;nlocal=0
  do z=anchor(3),last(3),2;do y=anchor(2),last(2),2;do x=anchor(1),last(1),2
    first=[x,y,z]
    if(.not.cg8_cell_eligible(first,s%origin,[c%nq,c%nr,c%ns],s%exclude_lower,s%exclude_upper)) cycle
    nlocal=nlocal+1;ijk=(first-anchor)/2+1;map(ijk(1),ijk(2),ijk(3))=nlocal
  enddo;enddo;enddo
  allocate(s%first(3,nlocal),s%owner(nlocal),s%owned_index(nlocal));nowned=0
  do z=anchor(3),last(3),2;do y=anchor(2),last(2),2;do x=anchor(1),last(1),2
    first=[x,y,z];ijk=(first-anchor)/2+1;j=map(ijk(1),ijk(2),ijk(3));if(j==0) cycle
    s%first(:,j)=first;r=node_rank(first,ranges,c%rank);s%owner(j)=r;s%owned_index(j)=0
    if(r==c%rank) then
      nowned=nowned+1;s%owned_index(j)=nowned
    endif
  enddo;enddo;enddo
  allocate(s%index(c%mq:c%pq,c%mr:c%pr,c%ms:c%ps));nfull=0;k=0
  do z=c%ms,c%ps;do y=c%mr,c%pr;do x=c%mq,c%pq
    point=[x,y,z];first=point-modulo(modulo(point,2)-modulo(s%origin,2),2)
    ijk=(first-anchor)/2+1;j=map(ijk(1),ijk(2),ijk(3))
    if(j==0) then
      nfull=nfull+1;s%index(x,y,z)=-nfull
    else
      s%index(x,y,z)=j;k=k+1
    endif
  enddo;enddo;enddo
  local_counts=[int(nowned,8),int(nfull,8),int(k,8)]
  call MPI_Allreduce(local_counts,global_counts,3,MPI_INTEGER8,MPI_SUM,c%comm,ierr);call mpi_check(ierr)
  if(global_counts(1)==0) call cgt_fail('Domain/buffers leave no complete projected interior cell')
  if(8_8*global_counts(1)/=global_counts(3)) call cgt_fail('Cell ownership does not cover exactly eight nodes')
  allocate(s%eta_cell(6,8,nowned),s%deta_cell(6,8,nowned),s%eta_full(6,8,nfull),s%deta_full(6,8,nfull))
  allocate(s%strain(6,8,nlocal),s%feedback(6,nlocal),s%seen(8,nlocal),s%full_seen(nfull))
  s%eta_cell=0;s%deta_cell=0;s%eta_full=0;s%deta_full=0
  s%strain=0;s%feedback=0;s%seen=.false.;s%full_seen=.false.
  counts=0
  do j=1,nlocal
    touched=.false.;r=s%owner(j)
    do k=1,8
      p=node_rank(slot_point(s%first(:,j),k),ranges,c%rank)
      if(r==c%rank.and.p/=c%rank) then
        counts(2,p)=counts(2,p)+1;touched(p)=.true.
      else if(r/=c%rank.and.p==c%rank) then
        counts(1,r)=counts(1,r)+1
      endif
    enddo
    if(r==c%rank) then
      where(touched) counts(3,:)=counts(3,:)+1
    else
      counts(4,r)=counts(4,r)+1
    endif
  enddo
  allocate(all_peers(0:np-1));positions=0
  do r=0,np-1
    all_peers(r)%rank=r;a=counts(1,r);b=counts(2,r);d=counts(3,r);e=counts(4,r)
    allocate(all_peers(r)%send_cell(a),all_peers(r)%send_slot(a),all_peers(r)%recv_cell(b),all_peers(r)%recv_slot(b))
    allocate(all_peers(r)%stress_send_cell(d),all_peers(r)%stress_recv_cell(e))
    allocate(all_peers(r)%strain_send(6,a),all_peers(r)%strain_recv(6,b))
    allocate(all_peers(r)%stress_send(6,d),all_peers(r)%stress_recv(6,e))
  enddo
  do j=1,nlocal
    touched=.false.;r=s%owner(j)
    do k=1,8
      p=node_rank(slot_point(s%first(:,j),k),ranges,c%rank)
      if(r==c%rank.and.p/=c%rank) then
        positions(2,p)=positions(2,p)+1;a=positions(2,p)
        all_peers(p)%recv_cell(a)=j;all_peers(p)%recv_slot(a)=k;touched(p)=.true.
      else if(r/=c%rank.and.p==c%rank) then
        positions(1,r)=positions(1,r)+1;a=positions(1,r)
        all_peers(r)%send_cell(a)=j;all_peers(r)%send_slot(a)=k
      endif
    enddo
    if(r==c%rank) then
      do p=0,np-1
        if(.not.touched(p)) cycle
        positions(3,p)=positions(3,p)+1;all_peers(p)%stress_send_cell(positions(3,p))=j
      enddo
    else
      positions(4,r)=positions(4,r)+1;all_peers(r)%stress_recv_cell(positions(4,r))=j
    endif
  enddo
  npeer=count(any(counts>0,dim=1));allocate(s%peers(npeer));a=0
  do r=0,np-1
    if(.not.any(counts(:,r)>0)) cycle
    a=a+1;s%peers(a)=all_peers(r)
  enddo
  ! Lists are ordered by global cell anchor then slot, independent of local ranges.
  call verify_routes(s)
  remote_slots=sum(counts(1,:))
  call MPI_Allreduce(remote_slots,global_remote,1,MPI_INTEGER,MPI_SUM,c%comm,ierr);call mpi_check(ierr)
  call MPI_Allreduce(npeer,max_peers,1,MPI_INTEGER,MPI_MAX,c%comm,ierr);call mpi_check(ierr)
  if(c%rank==0) then
    write(*,'(A,I0,A,I0)') 'CG-T remote strain slots=',global_remote,', max cell peers=',max_peers
    write(*,'(A,3I12)') 'CG-T owned cells/full nodes/projected nodes=',global_counts
    write(*,'(A,I0)') 'CG-T memory-variable bytes=',768_8*(global_counts(1)+global_counts(2))
  endif
 end subroutine
 subroutine verify_routes(s)
  type(cgt_state),intent(in)::s
  integer::p,i,j,k,n,ierr
  integer,allocatable::send(:,:),recv(:,:)
  do p=1,size(s%peers)
    allocate(send(4,size(s%peers(p)%send_cell)),recv(4,size(s%peers(p)%recv_cell)))
    do i=1,size(send,2)
      j=s%peers(p)%send_cell(i);send(:,i)=[s%first(:,j),s%peers(p)%send_slot(i)]
    enddo
    call MPI_Sendrecv(send,size(send),MPI_INTEGER,s%peers(p)%rank,20041, &
      recv,size(recv),MPI_INTEGER,s%peers(p)%rank,20041,s%comm,MPI_STATUS_IGNORE,ierr);call mpi_check(ierr)
    do i=1,size(recv,2)
      j=s%peers(p)%recv_cell(i);k=s%peers(p)%recv_slot(i)
      if(any(recv(:,i)/=[s%first(:,j),k])) call cgt_fail('Nonmatching global cell/slot strain routes')
    enddo
    deallocate(send,recv)
    allocate(send(3,size(s%peers(p)%stress_send_cell)),recv(3,size(s%peers(p)%stress_recv_cell)))
    do i=1,size(send,2)
      send(:,i)=s%first(:,s%peers(p)%stress_send_cell(i))
    enddo
    call MPI_Sendrecv(send,size(send),MPI_INTEGER,s%peers(p)%rank,20042, &
      recv,size(recv),MPI_INTEGER,s%peers(p)%rank,20042,s%comm,MPI_STATUS_IGNORE,ierr);call mpi_check(ierr)
    do i=1,size(recv,2)
      if(any(recv(:,i)/=s%first(:,s%peers(p)%stress_recv_cell(i)))) call cgt_fail('Nonmatching memory stress routes')
    enddo
    deallocate(send,recv)
  enddo
 end subroutine
 subroutine exchange_cgt_stage(s)
  type(cgt_state),intent(inout)::s
  integer::p,i,j,k,n,ierr,nr
  integer,allocatable::requests(:)
  allocate(requests(4*size(s%peers)));nr=0
  do j=1,size(s%owner)
    i=s%owned_index(j);if(i==0) cycle
    s%feedback(:,j)=0
    do k=1,8
      s%feedback(:,j)=s%feedback(:,j)+s%eta_cell(:,k,i)
    enddo
  enddo
  do p=1,size(s%peers)
    n=size(s%peers(p)%strain_recv)
    if(n>0) then
      nr=nr+1
      call MPI_Irecv(s%peers(p)%strain_recv,n,MPI_DOUBLE_PRECISION,s%peers(p)%rank,20043, &
        s%comm,requests(nr),ierr);call mpi_check(ierr)
    endif
    n=size(s%peers(p)%stress_recv)
    if(n>0) then
      nr=nr+1
      call MPI_Irecv(s%peers(p)%stress_recv,n,MPI_DOUBLE_PRECISION,s%peers(p)%rank,20044, &
        s%comm,requests(nr),ierr);call mpi_check(ierr)
    endif
  enddo
  do p=1,size(s%peers)
    do i=1,size(s%peers(p)%send_cell)
      j=s%peers(p)%send_cell(i);k=s%peers(p)%send_slot(i)
      if(.not.s%seen(k,j)) call cgt_fail('Missing outgoing strain slot')
      s%peers(p)%strain_send(:,i)=s%strain(:,k,j)
    enddo
    do i=1,size(s%peers(p)%stress_send_cell)
      j=s%peers(p)%stress_send_cell(i);s%peers(p)%stress_send(:,i)=s%feedback(:,j)
    enddo
    n=size(s%peers(p)%strain_send)
    if(n>0) then
      nr=nr+1
      call MPI_Isend(s%peers(p)%strain_send,n,MPI_DOUBLE_PRECISION,s%peers(p)%rank,20043, &
        s%comm,requests(nr),ierr);call mpi_check(ierr)
    endif
    n=size(s%peers(p)%stress_send)
    if(n>0) then
      nr=nr+1
      call MPI_Isend(s%peers(p)%stress_send,n,MPI_DOUBLE_PRECISION,s%peers(p)%rank,20044, &
        s%comm,requests(nr),ierr);call mpi_check(ierr)
    endif
  enddo
  if(nr>0) then
    call MPI_Waitall(nr,requests,MPI_STATUSES_IGNORE,ierr);call mpi_check(ierr)
  endif
  do p=1,size(s%peers)
    do i=1,size(s%peers(p)%recv_cell)
      j=s%peers(p)%recv_cell(i);k=s%peers(p)%recv_slot(i)
      if(s%seen(k,j)) call cgt_fail('Duplicate remote strain slot')
      s%strain(:,k,j)=s%peers(p)%strain_recv(:,i);s%seen(k,j)=.true.
    enddo
    do i=1,size(s%peers(p)%stress_recv_cell)
      j=s%peers(p)%stress_recv_cell(i);s%feedback(:,j)=s%peers(p)%stress_recv(:,i)
    enddo
  enddo
 end subroutine
end module
