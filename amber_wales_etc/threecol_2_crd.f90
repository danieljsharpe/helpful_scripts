! Fortran script to read a coordinates file in the 3-column plain x-y-z format dumped by EXTRACTMIN
! and EXTRACTTS keywords of PATHSAMPLE and write the coordinates in AMBER .crd format

! compile with: gfortran threecol_2_crd.f90 -o threecol_2_crd

! Daniel J. Sharpe, Mar 2019

program threecol_2_crd

integer :: natoms, i
character(len=32) :: arg
double precision, allocatable :: coords(:,:)

CALL get_command_argument(1,arg)
read(arg,"(I4)") natoms

allocate(coords(3,natoms))

open(12,file="extractedmin")

read(unit=12,fmt="(3F25.15)") coords

!do i=1,3
!    print *, coords(1,i), coords(2,i), coords(3,i)
!end do

write(*,fmt="(10F8.3)") (coords(1,i),coords(2,i),coords(3,i),i=1,natoms)

close(12)

deallocate(coords)

end program threecol_2_crd
