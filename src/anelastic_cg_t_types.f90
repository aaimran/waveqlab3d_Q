module anelastic_cg_t_types
 use common,only:wp
 use anelastic_cg8_model,only:cg8_coefficients
 implicit none
 private
 public::cgt_state,cgt_peer
 type cgt_peer
   integer::rank=0
   integer,allocatable::send_cell(:),send_slot(:),recv_cell(:),recv_slot(:)
   integer,allocatable::stress_send_cell(:),stress_recv_cell(:)
   real(wp),allocatable::strain_send(:,:),strain_recv(:,:),stress_send(:,:),stress_recv(:,:)
 end type
 type cgt_state
   type(cg8_coefficients)::coeff
   integer::comm=0,rank=0,origin(3)=1,exclude_lower(3)=0,exclude_upper(3)=0
   integer,allocatable::index(:,:,:),first(:,:),owner(:),owned_index(:)
   type(cgt_peer),allocatable::peers(:)
   real(wp),allocatable::eta_cell(:,:,:),deta_cell(:,:,:),eta_full(:,:,:),deta_full(:,:,:)
   real(wp),allocatable::strain(:,:,:),feedback(:,:)
   logical,allocatable::seen(:,:),full_seen(:)
 end type
end module
