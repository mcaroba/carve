!----------------------------------------------------------------------!
!
! carve.f90 (this code) is Copyright (c) 2013-2020 of Miguel A. Caro
! and was written during his employment at Aalto University.
!
! The code is released under the Creative Commons Attribution
! Share-Alike International License:
!
! https://creativecommons.org/licenses/by-sa/4.0/
!
! If you want to use the code under different terms than those covered
! by the license, or have other questions, bug reports, etc., you
! can contact Miguel Caro:
!
! miguel.caro@aalto.fi
! mcaroba@gmail.com
!
! The official repository of this code is on Github:
!
! https://github.com/mcaroba/carve
!
! The proper way of giving attribution when using the code is to cite
! the following papers (citation info will be updated when a better
! reference is available):
!
! M.A. Caro, R. Zoubkoff, O. Lopez-Acevedo and T. Laurila
! "Atomic and electronic structure of tetrahedral amorphous carbon
! surfaces from density functional theory: Properties and simulation
! strategies"
! Carbon 77, 1168 (2014)
! https://doi.org/10.1016/j.carbon.2014.06.060
!
! M.A. Caro, G. Csányi, T. Laurila and V.L. Deringer
! "Machine learning driven simulated deposition of carbon films: From
! low-density to diamondlike amorphous carbon"
! Phys. Rev. B 102, 174201 (2020)
! https://doi.org/10.1103/PhysRevB.102.174201.
!----------------------------------------------------------------------!






!----------------------------------------------------------------------!
module subroutines
! Several subroutines

  implicit none

  contains

  function cross(u,v)
! Obtains the cross product of vectors u and v
    real*8 :: cross(1:3)
    real*8, intent(in) :: u(1:3), v(1:3)

    cross = (/ &
      u(2)*v(3)-u(3)*v(2), &
      u(3)*v(1)-u(1)*v(3), &
      u(1)*v(2)-u(2)*v(1) &
      /)
  end function cross
!----------------------------------------------------------------------!
  function norm(u)
! Obtains the norm of u
    real*8 :: norm
    real*8, intent(in) :: u(1:3)

    norm = dsqrt(u(1)**2+u(2)**2+u(3)**2)
  end function norm
!----------------------------------------------------------------------!
  function dot(u,v)
! Obtains the dot product of vectors u and v
    real*8 :: dot
    real*8, intent(in) :: u(1:3), v(1:3)

    dot = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
  end function dot
!----------------------------------------------------------------------!
! This subroutine returns the distance between ri and rj under certain boundary conditions
  subroutine get_distance(posi, posj, a, b, c, PBC, dist, d)

    implicit none

    real*8, intent(in) :: posi(1:3), posj(1:3), a(1:3), b(1:3), c(1:3)
    logical, intent(in) :: PBC(1:3)
    real*8, intent(out) :: d
    real*8, intent(out) :: dist(1:3)
    real*8 :: d2, L(1:3)
    real*8 :: mat(1:3,1:3), md, indices_real(1:3), res, res_opt, dist_opt(1:3), dist_temp(1:3)
    real*8, save :: a0(1:3) = 0.d0, b0(1:3) = 0.d0, c0(1:3) = 0.d0, mat_inv(1:3,1:3) = 0.d0
    integer :: i, j, k, indices(1:3), i_opt, j_opt, k_opt
    logical :: lattice_check_a(1:3), lattice_check_b(1:3), lattice_check_c(1:3)

    if( a(2) == 0.d0 .and. a(3) == 0.d0 .and. b(1) == 0.d0 .and. &
        b(3) == 0.d0 .and. c(1) == 0.d0 .and. c(2) == 0.d0 )then
!     Fast solution for orthorhombic cells
      L = (/ a(1), b(2), c(3) /)
      d2 = 0.d0
      do i = 1, 3
        if( PBC(i) )then
          dist(i) = modulo(posj(i) - posi(i), L(i))
          if( dist(i) > L(i)/2.d0 )then
            dist(i) = dist(i) - L(i)
          end if
        else
          dist(i) = posj(i) - posi(i)
        end if
        d2 = d2 + dist(i)**2
      end do
      d = dsqrt(d2)
    else if( all( PBC ) )then
!     Slow solution for other unit cells
      lattice_check_a = ( a /= a0 )
      lattice_check_b = ( b /= b0 )
      lattice_check_c = ( c /= c0 )
      if( any(lattice_check_a) .or. any(lattice_check_b) .or. any(lattice_check_c) )then
        a0 = a
        b0 = b
        c0 = c
!       We construct our matrix to get the MIC only if the lattice vectors have changed
        mat(1,1) = dot_product(a, a)
        mat(1,2) = dot_product(a, b)
        mat(1,3) = dot_product(a, c)
        mat(2,1) = mat(1,2)
        mat(2,2) = dot_product(b, b)
        mat(2,3) = dot_product(b, c)
        mat(3,1) = mat(1,3)
        mat(3,2) = mat(2,3)
        mat(3,3) = dot_product(c, c)
!       We compute the inverse of this matrix analytically
        md = -mat(1,3)**2*mat(2,2) + 2.d0*mat(1,2)*mat(1,3)*mat(2,3) - mat(1,1)*mat(2,3)**2 &
             - mat(1,2)**2*mat(3,3) + mat(1,1)*mat(2,2)*mat(3,3)
        mat_inv(1,1) = mat(2,2)*mat(3,3) - mat(2,3)**2
        mat_inv(1,2) = mat(1,3)*mat(2,3) - mat(1,2)*mat(3,3)
        mat_inv(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
        mat_inv(2,1) = mat_inv(1,2)
        mat_inv(2,2) = mat(1,1)*mat(3,3) - mat(1,3)**2
        mat_inv(2,3) = mat(1,2)*mat(1,3) - mat(1,1)*mat(2,3)
        mat_inv(3,1) = mat_inv(1,3)
        mat_inv(3,2) = mat_inv(2,3)
        mat_inv(3,3) = mat(1,1)*mat(2,2) - mat(1,2)**2
        mat_inv = mat_inv / md
      end if
      dist = posj - posi
      indices_real = -2.d0 * (/ dot_product(dist, a), dot_product(dist, b), dot_product(dist,c) /)
      indices_real = matmul(mat_inv, indices_real)
!     Closest integer solution
      indices(1:3) = nint(indices_real(1:3))
!     We bruteforce the integer solution among the 27 points surrounding the real solution
      res_opt = 1.d10
      do i = indices(1)-1, indices(1)+1
        do j = indices(2)-1, indices(2)+1
          do k = indices(3)-1, indices(3)+1
            dist_temp(1:3) = dist(1:3) + dfloat(i)*a(1:3) + dfloat(j)*b(1:3) + dfloat(k)*c(1:3)
            res = dot_product(dist_temp, dist_temp)
            if( res < res_opt )then
              res_opt = res
              dist_opt = dist_temp
            end if
          end do
        end do
      end do
      dist = dist_opt
      d = dsqrt( dot_product(dist,dist) )
    else
      write(*,*) "Sorry, non-orthorhombic unit cells only work in combination with full PBC"
    end if

  return
  end subroutine get_distance
!----------------------------------------------------------------------!
end module subroutines















!----------------------------------------------------------------------!
!                   MAIN PROGRAM STARTS HERE                           !
!----------------------------------------------------------------------!
program carve

  use subroutines

! This program processes CONTCAR files in VASP and
! carves out a chunk off the system. It should probably be totally rewritten
! since it's fairly shitty at the moment, but it does the job

! ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
!
! ATTENTION: This code only works for orthorhombic cells <- I'm actually not sure if this is true
!
! ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION

  implicit none

  integer :: error, counter, i, j, types, n_atoms, n_points, k, i2, j2, k2, ia, ib, ic
  integer :: i3, j3, ci(4,3), i4
  character*16384 :: file, elements, header
  real*8 :: sigma, cutoff, f, a(3), b(3), c(3), axb(3), S, temp(3), dz, z, mass, convfactor
  real*8, allocatable :: atom_position(:,:), atom_mass(:), density(:), smeared_density(:)
  real*8, allocatable :: position_array(:), density_sp3(:), density_sp2(:)
  logical :: smearing = .TRUE.
  character*2, allocatable :: atom_species(:), atom_types(:)
  integer, allocatable :: atom_number(:), atom_coordination(:), atom_neighbors(:,:)
  character*2 :: cjunk
  character*1 :: mode
  real*8 :: nn_cutoff, mass_sp2, mass_sp3, expectation, occurrence, thickness = 3.d0
  integer :: N_sp2, N_sp3, N_chains4
  integer :: l, npoints, natoms_in_slice
  integer, allocatable :: qi(:)
  real*8, allocatable :: qz(:,:)
  real*8 :: v1(3), v2(3), v3(3), pi, z_step, sum1, sum2, zmax, vector(3)
  real*8 :: dHC = 0.9d0
  real*8 :: pos(1:3), d, dmin
  integer :: imin, jmin, kmin
  integer :: center_atom = 1, ia2, ib2, ic2
  real*8 :: r, r2
  real*8 :: dHH_min = 1.1d0, pos2(1:3)
  real*8, allocatable :: new_H_pos(:,:), new_H_pos_temp(:,:)
  integer :: n_new_H, n_keep
  logical :: too_close, keep_running, fix_atoms, override_ortho = .false.
  logical, allocatable :: remove_atom(:)
  real*8 :: dist(1:3), dist2(1:3), fix_in, fix_out, shift(1:3)
  character*64 :: output_format = "xyz", argument
! This is for VASP output. At the moment we handle structures with C, O and H, in that order
  integer :: n_vasp(1:3) = 0
  real*8 :: x_max, x_min, y_max, y_min, z_max, z_min, nn_cutoff_max = 1.9d0

  pi = dacos(-1.d0)


! This is now given as a function at the end of the file
!  nn_cutoff = 1.9d0


!  read(*,*, iostat=error) thickness, center_atom, dHC, dHH_min, file, output_format, fix_in

  do i = 1, 8
    call getarg(i, argument)
    if( argument(1:5) == "rcut=" )then
      read(argument(6:64), *) thickness
    else if( argument(1:13) == "central_atom=" )then
      read(argument(14:64), *) center_atom
    else if( argument(1:4) == "dHC=" )then
      read(argument(5:64), *) dHC
    else if( argument(1:8) == "dHH_min=" )then
      read(argument(9:64), *) dHH_min
    else if( argument(1:11) == "atoms_file=" )then
      read(argument(12:64), *) file
    else if( argument(1:7) == "format=" )then
      read(argument(8:64), *) output_format
    else if( argument(1:5) == "rfix=" )then
      read(argument(6:64), *) fix_in
    else if( argument(1:5) == "override_ortho" )then
      override_ortho = .true.
    else if( argument == "" )then
      continue
    else
      write(0,*) "ERROR: I don't understand command line argument ", argument
      stop
    end if
  end do

  shift(1:3) = (/ 3.d0+thickness, 3.d0+thickness, 3.d0+thickness /)

!  output_format = trim(adjustl( output_format ))

  if( output_format == "fix_in" )then
    fix_out = fix_in
    fix_in = 0.d0
    output_format = "vasp"
  else if( output_format == "fix_out" )then
    fix_out = 1.d10
    fix_in = fix_in
    output_format = "vasp"
  else
    fix_out = 1.d10
    fix_in = 1.d10
  end if


  open(unit=10, file=file, status='old')

! Read CONTCAR
!
! Header
  read(10, '(A)', iostat=error) header
!
! Prefactor
  read(10,*) f
!
! Lattice Vectors
  read(10,*) a(1), a(2), a(3)
  read(10,*) b(1), b(2), b(3)
  read(10,*) c(1), c(2), c(3)
  a = f*a
  b = f*b
  c = f*c

  if( ( (a(2) /= 0.d0) .or. (a(3) /= 0.d0) .or. (b(1) /= 0.d0) .or. (b(3) /= 0.d0) .or. &
      (c(1) /= 0.d0) .or. (c(2) /= 0.d0) ) .and. ( .not. override_ortho ) )then
    write(0, *) "ERROR: only orthorhombic unit cells can currently be processed!"
    stop
  end if

!
! Check atomic types present
  read(10,'(A)') elements
! If the elements line contains numbers, then an old-format POSCAR is being used. In
! such case, we assume the elements are given in the header, which should in principle
! be a comment line but is often used for that purpose
  read(elements, *, iostat=error) i
  if( error == 0 )then
    elements = header
    backspace(10)
  end if
  counter=0
  error=0
  do while ( error == 0)
    read(elements,*,iostat=error) ( cjunk, i=0,counter )
    counter=counter+1
  end do
  types = counter - 1
  allocate ( atom_types(1:types) )
  read(elements, *) ( atom_types(i), i=1,types )
!
! Read number of atoms of each type
  allocate ( atom_number(1:types) )
  read(10,*) ( atom_number(i), i=1,types )
  n_atoms = 0
  do i=1,types
    n_atoms = n_atoms + atom_number(i)
  end do
!
! Read in atomic positions
  allocate ( atom_species(1:n_atoms) )
  allocate ( atom_position(1:n_atoms,1:3) )
  counter = 1
  do i=1,types
    do j=1,atom_number(i)
      atom_species(counter) = trim(adjustl( atom_types(i) ))
      counter = counter + 1
    end do
  end do
  read(10,*) mode
  if(mode == 's' .or. mode == 'S')then
    read(10,*) mode
  end if
  do i=1,n_atoms
    read(10,*) ( atom_position(i,j), j=1,3 )
  end do
  if( mode == 'd' .or. mode == 'D')then
    do i=1,n_atoms
      temp(1) = a(1)*atom_position(i,1) + b(1)*atom_position(i,2) + c(1)*atom_position(i,3)
      temp(2) = a(2)*atom_position(i,1) + b(2)*atom_position(i,2) + c(2)*atom_position(i,3)
      temp(3) = a(3)*atom_position(i,1) + b(3)*atom_position(i,2) + c(3)*atom_position(i,3)
      do j=1,3
        atom_position(i,j) = temp(j)
      end do
    end do
  end if
!
  close(10)
!
! Obtain coordination
!
! This is an old and inefficient brute-force code
!
  allocate ( atom_coordination(1:n_atoms) )
  allocate ( atom_neighbors(1:n_atoms, 1:10) )
  atom_neighbors = 0
  do i=1,n_atoms
    counter=0
    do j=1,n_atoms
      call get_distance(atom_position(j, 1:3), atom_position(i, 1:3), a, b, c, (/ .true., .true., .true. /), dist, d)
      if( i /= j .and. d < nn_cutoff(atom_species(i), atom_species(j)) )then
        counter = counter + 1
        atom_neighbors(i,counter) = j
      end if
    end do
    atom_coordination(i) = counter
  end do


! Determine which atoms are in the slice
!
! Check for each atom in the slice if and how many of its neighbors are outside the slice
! For each neighbor outside the slice place a H atom dHC A from the atoms in the slice
! along the vector connecting it with the neighbor outside
!
  allocate( remove_atom(1:n_atoms) )
  remove_atom = .true.
  keep_running = .true.
  do while( keep_running )
    keep_running = .false.
    do i=1,n_atoms
      call get_distance(atom_position(center_atom, 1:3), atom_position(i, 1:3), a, b, c, &
                        (/ .true., .true., .true. /), dist, d)
      if( .not. remove_atom(i) )then
        cycle
      else
        if( d <= thickness .and. atom_species(i) == "C" )then
          remove_atom(i) = .false.
          keep_running = .true.
        else if( atom_species(i) == "H" .or. atom_species(i) == "O" )then
          do j = 1, atom_coordination(i)
            k = atom_neighbors(i, j)
            if( .not. remove_atom(k) )then
              remove_atom(i) = .false.
              keep_running = .true.
              do j2 = 1, atom_coordination(i)
                k2 = atom_neighbors(i, j2)
                remove_atom(k2) = .false.
              end do
              cycle
            end if
          end do
        end if
      end if
    end do
  end do


! Check that the unit cell is big enough for the procedure to work. Else,
! produce an error message
  x_max = -1.d10
  x_min = 1.d10
  y_max = -1.d10
  y_min = 1.d10
  z_max = -1.d10
  z_min = 1.d10
  do i = 1, n_atoms
    if( .not. remove_atom(i) )then
      call get_distance(atom_position(center_atom, 1:3), atom_position(i, 1:3), a, b, c, &
                        (/ .true., .true., .true. /), dist, d)
      if( dist(1) > x_max )then
        x_max = dist(1)
      end if
      if( dist(2) > y_max )then
        y_max = dist(2)
      end if
      if( dist(3) > z_max )then
        z_max = dist(3)
      end if
      if( dist(1) < x_min )then
        x_min = dist(1)
      end if
      if( dist(2) < y_min )then
        y_min = dist(2)
      end if
      if( dist(3) < z_min )then
        z_min = dist(3)
      end if
    end if
  end do
  if( x_max - x_min >= a(1) - nn_cutoff_max )then
    write(0, *) "WARNING: your unit cell may be too small along the x direction!"
    write(0, *) "Make sure the results make sense (e.g., keep an eye for radicals)."
    write(0, *) "You can use ASE to easily make it bigger, e.g.,"
    write(0, *) "atoms *= (2,1,1)"
  end if
  if( y_max - y_min >= b(2) - nn_cutoff_max )then
    write(0, *) "WARNING: your unit cell may be too small along the y direction!"
    write(0, *) "Make sure the results make sense (e.g., keep an eye for radicals)."
    write(0, *) "You can use ASE to easily make it bigger, e.g.,"
    write(0, *) "atoms *= (1,2,1)"
  end if
  if( z_max - z_min >= c(3) - nn_cutoff_max )then
    write(0, *) "WARNING: your unit cell may be too small along the z direction!"
    write(0, *) "Make sure the results make sense (e.g., keep an eye for radicals)."
    write(0, *) "You can use ASE to easily make it bigger, e.g.,"
    write(0, *) "atoms *= (1,1,2)"
  end if


! Figure out how many atoms we keep
  n_keep = 0
  n_new_H = 0
  do i=1,n_atoms
    if( .not. remove_atom(i) )then
      n_keep = n_keep + 1
    end if
    if( atom_species(i) == "C" .and. (.not. remove_atom(i)) )then
      do j = 1, atom_coordination(i)
        k = atom_neighbors(i, j)
        if( atom_species(k) == "C" .and. remove_atom(k) )then
          n_new_H = n_new_H + 1
        end if
      end do
    end if
  end do
  allocate( new_H_pos(1:6,1:n_new_H) )

  if( output_format == "xyz" )then
    write(*,*) n_keep + n_new_H
    write(*,*) "comment=""This passivated cluster structure was generated with Miguel Caro's carve.f90 code"""
  end if

! Now store info to set up the passivating Hs
  n_new_H = 0
! We print the central atom (where the core is located) first for convenience
  call get_distance(atom_position(center_atom, 1:3), atom_position(center_atom, 1:3), a, b, c, &
                    (/ .true., .true., .true. /), dist, d)
! I'm sure this bit of code can be written much better
  if( output_format == "xyz" )then
    write(*,*) atom_species(center_atom), dist(:)
  else if( output_format == "vasp" )then
    if( atom_species(center_atom) == "C" )then
      n_vasp(1) = n_vasp(1) + 1
    else if( atom_species(center_atom) == "O" )then
      n_vasp(2) = n_vasp(2) + 1
    else if( atom_species(center_atom) == "H" )then
      n_vasp(3) = n_vasp(3) + 1
    end if
  end if
  do i=1,n_atoms
    if( (.not. remove_atom(i)) .and. (i /= center_atom) )then
      call get_distance(atom_position(center_atom, 1:3), atom_position(i, 1:3), a, b, c, &
                        (/ .true., .true., .true. /), dist, d)
      if( output_format == "xyz" )then
        write(*,*) atom_species(i), dist(:)
      else if( output_format == "vasp" )then
        if( atom_species(i) == "C" )then
          n_vasp(1) = n_vasp(1) + 1
        else if( atom_species(i) == "O" )then
          n_vasp(2) = n_vasp(2) + 1
        else if( atom_species(i) == "H" )then
          n_vasp(3) = n_vasp(3) + 1
        end if
      end if
    end if
    if( atom_species(i) == "C" .and. (.not. remove_atom(i)) )then
      do j = 1, atom_coordination(i)
        k = atom_neighbors(i, j)
        if( atom_species(k) == "C" .and. remove_atom(k) )then
          n_new_H = n_new_H + 1
          call get_distance(atom_position(center_atom, 1:3), atom_position(i, 1:3), a, b, c, &
                            (/ .true., .true., .true. /), dist, d)
          new_H_pos(1:3,n_new_H) = dist(1:3)
          call get_distance(atom_position(i, 1:3), atom_position(k, 1:3), a, b, c, &
                            (/ .true., .true., .true. /), dist, d)
          new_H_pos(4:6,n_new_H) = dist(1:3)
        end if
      end do
    end if
  end do



  if( output_format == "vasp" )then
    n_vasp(3) = n_vasp(3) + n_new_H
    header = ""
    elements = ""
    if( n_vasp(1) > 0 )then
      write(header, '(A,1X,A)') trim(adjustl(header)), "C"
      write(elements, '(A,1X,I6)') trim(adjustl(elements)), n_vasp(1)
    end if
    if( n_vasp(2) > 0 )then
      write(header, '(A,1X,A)') trim(adjustl(header)), "O"
      write(elements, '(A,1X,I6)') trim(adjustl(elements)), n_vasp(2)
    end if
    if( n_vasp(3) > 0 )then
      write(header, '(A,1X,A)') trim(adjustl(header)), "H"
      write(elements, '(A,1X,I6)') trim(adjustl(elements)), n_vasp(3)
    end if
    write(*,*) "This passivated cluster structure was generated with Miguel Caro's carve.f90 code"
    write(*,*) 1.d0
    write(*,*) 2.d0*shift(1), 0.d0, 0.d0
    write(*,*) 0.d0, 2.d0*shift(2), 0.d0
    write(*,*) 0.d0, 0.d0, 2.d0*shift(3)
    write(*,*) trim(adjustl(header))
    write(*,*) trim(adjustl(elements))
    write(*,'(A)') "Selective Dynamics"
    write(*,'(A)') "Cartesian"
    do i=1,n_atoms
      if( (.not. remove_atom(i)) .and. (atom_species(i) == "C") )then
        call get_distance(atom_position(center_atom, 1:3), atom_position(i, 1:3), a, b, c, &
                          (/ .true., .true., .true. /), dist, d)
        if( d >= fix_in .and. d < fix_out )then
          write(*, *) dist(1:3)+shift(1:3), "  F   F   F"
        else
          write(*, *) dist(1:3)+shift(1:3), "  T   T   T"
        end if
      end if
    end do
    do i=1,n_atoms
      if( (.not. remove_atom(i)) .and. (atom_species(i) == "O") )then
        call get_distance(atom_position(center_atom, 1:3), atom_position(i, 1:3), a, b, c, &
                          (/ .true., .true., .true. /), dist, d)
        if( d >= fix_in .and. d < fix_out )then
          write(*, *) dist(1:3)+shift(1:3), "  F   F   F"
        else
          write(*, *) dist(1:3)+shift(1:3), "  T   T   T"
        end if
      end if
    end do
    do i=1,n_atoms
      if( (.not. remove_atom(i)) .and. (atom_species(i) == "H") )then
        call get_distance(atom_position(center_atom, 1:3), atom_position(i, 1:3), a, b, c, &
                          (/ .true., .true., .true. /), dist, d)
        if( d >= fix_in .and. d < fix_out )then
          write(*, *) dist(1:3)+shift(1:3), "  F   F   F"
        else
          write(*, *) dist(1:3)+shift(1:3), "  T   T   T"
        end if
      end if
    end do
  end if






! Check that H atoms are not too close to each other
  too_close = .true.
  do while(too_close)
    too_close = .false.
    do i = 1, n_new_H
      do j = i+1, n_new_H
        pos(1:3) = new_H_pos(1:3,i) + dHC * new_H_pos(4:6,i) / norm(new_H_pos(4:6,i))
        pos2(1:3) = new_H_pos(1:3,j) + dHC * new_H_pos(4:6,j) / norm(new_H_pos(4:6,j))
        d = dsqrt( (pos(1) - pos2(1))**2 + &
                   (pos(2) - pos2(2))**2 + &
                   (pos(3) - pos2(3))**2 )
        if( d < dHH_min )then
          new_H_pos(4:6,i) = new_H_pos(4:6,i) + 0.05d0 * (pos(1:3) - pos2(1:3))
          new_H_pos(4:6,j) = new_H_pos(4:6,j) + 0.05d0 * (pos2(1:3) - pos(1:3))
          too_close = .true.
        end if
      end do
    end do
  end do


  do i = 1, n_new_H
    pos(1:3) = new_H_pos(1:3,i) + dHC * new_H_pos(4:6,i) / norm(new_H_pos(4:6,i))
    if( output_format == "xyz" )then
      write(*, *) "H", pos(1:3)
    else if( output_format == "vasp" )then
      d = dsqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
      if( d >= fix_in .and. d < fix_out )then
        write(*, *) pos(1:3)+shift(1:3), "  F   F   F"
      else
        write(*, *) pos(1:3)+shift(1:3), "  T   T   T"
      end if
    end if
  end do


end program carve







function nn_cutoff(species1, species2)

  implicit none

  character*2 :: species1, species2
  real*8 :: nn_cutoff

  if( trim(adjustl(species1)) == "C" )then
    if( trim(adjustl(species2)) == "C" )then
      nn_cutoff = 1.9d0
      return
    else if( trim(adjustl(species2)) == "H" )then
      nn_cutoff = 1.3d0
      return
    else if( trim(adjustl(species2)) == "O" )then
      nn_cutoff = 1.7d0
      return
    end if
  else if( trim(adjustl(species1)) == "H" )then
    if( trim(adjustl(species2)) == "C" )then
      nn_cutoff = 1.3d0
      return
    else if( trim(adjustl(species2)) == "H" )then
      nn_cutoff = 1.d0
      return  
    else if( trim(adjustl(species2)) == "O" )then
      nn_cutoff = 1.3d0
      return
    end	if
  else if( trim(adjustl(species1)) == "O" )then
    if( trim(adjustl(species2)) == "C" )then
      nn_cutoff = 1.7d0
      return
    else if( trim(adjustl(species2)) == "H" )then
      nn_cutoff = 1.3d0
      return  
    else if( trim(adjustl(species2)) == "O" )then
      nn_cutoff = 1.5d0
      return
    end	if
  end if

  return

end function
