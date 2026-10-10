module anelastic_cg8_layout
 implicit none
 private
 public::cg8_mechanism,cg8_cell_eligible,cg8_minimum_buffer
contains
 pure integer function cg8_mechanism(point,origin) result(k)
  integer,intent(in)::point(3),origin(3)
  k=1+modulo(modulo(point(1),2)-modulo(origin(1),2),2)+ &
    2*modulo(modulo(point(2),2)-modulo(origin(2),2),2)+4*modulo(modulo(point(3),2)-modulo(origin(3),2),2)
 end function
 pure logical function cg8_cell_eligible(point,origin,shape,lower,upper) result(eligible)
  integer,intent(in)::point(3),origin(3),shape(3),lower(3),upper(3)
  integer::first(3)
  first=point-modulo(modulo(point,2)-modulo(origin,2),2)
  eligible=all(first>lower).and.all(first+1<=shape-upper)
 end function
 pure integer function cg8_minimum_buffer(closure,reach) result(n)
  integer,intent(in)::closure,reach
  n=2*((closure+reach+3)/2)
 end function
end module
