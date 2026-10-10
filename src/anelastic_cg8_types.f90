module anelastic_cg8_types
 use common,only:wp
 use anelastic_cg8_model,only:cg8_coefficients
 implicit none
 type :: cg8_state
   type(cg8_coefficients)::coarse,full
   integer,allocatable::index(:,:,:)
   real(wp),allocatable::eta_cg(:,:),deta_cg(:,:),eta_full(:,:,:),deta_full(:,:,:)
   integer::origin(3)=1,exclude_lower(3)=0,exclude_upper(3)=0
 end type
end module
