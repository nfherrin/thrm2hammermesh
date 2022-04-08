PROGRAM thrm2hammer
  IMPLICIT NONE
  !globals are up here as is the driver
  INTEGER(8) :: num_points,num_tets
  CHARACTER(64) :: mesh_file
  INTEGER(8),ALLOCATABLE :: tetgeom(:,:)
  REAL(8),ALLOCATABLE :: point(:,:),tetvol(:)
  REAL(8),PARAMETER :: PI=4.D0*DATAN(1.D0)
  INTEGER(8),PARAMETER::in_unit=20,out_unit=30

  CALL readin_cmdline()

  CALL readin_mesh()

  CALL computevol()

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE readin_cmdline()
    CHARACTER(64) :: temp_string
    INTEGER(8) :: arg_count

    !get arguments for the 3 potential inputs
    arg_count = COMMAND_ARGUMENT_COUNT()
    IF(arg_count .LT. 1)STOP 'no arguments!'
    IF(arg_count .GT. 1)STOP 'too many arguments! only need thrm file!'
    CALL GET_COMMAND_ARGUMENT(1, temp_string)
    mesh_file=temp_string

  ENDSUBROUTINE readin_cmdline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE readin_mesh()
    INTEGER(8) :: i,tempint
    !read in the mesh file
    OPEN(UNIT = in_unit, FILE = mesh_file, ACTION = "read", IOSTAT=tempint)
    IF(tempint .NE. 0)STOP 'error opening mesh file'

    READ(in_unit,*)num_points
    READ(in_unit,*)num_tets
    READ(in_unit,*)
    READ(in_unit,*)

    ALLOCATE(point(num_points,3),tetgeom(num_tets,4))

    DO i=1,num_points
      READ(in_unit,*)tempint,point(i,:)
      IF(i .NE. tempint)STOP 'points out of order'
    ENDDO

    DO i=1,num_tets
      READ(in_unit,*)tempint
      IF(i .NE. tempint)STOP 'element regions out of order'
    ENDDO

    DO i=1,num_tets
      READ(in_unit,*)tempint,tetgeom(i,:)
      IF(i .NE. tempint)STOP 'element geometry out of order'
    ENDDO

    CLOSE(in_unit)
  ENDSUBROUTINE readin_mesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE computevol()
    INTEGER(8) :: i
    REAL(8) :: totalvol1
    REAL(8) :: a(3),b(3),c(3),d(3)

    ALLOCATE(tetvol(num_tets))
    tetvol(:)=0
    totalvol1=0
    !compute the tet volumes
    DO i=1,num_tets
      a(:)=point(tetgeom(i,1),:)
      b(:)=point(tetgeom(i,2),:)
      c(:)=point(tetgeom(i,3),:)
      d(:)=point(tetgeom(i,4),:)
      tetvol(i)=ABS((-c(2)*d(1)+b(2)*(-c(1)+d(1))+b(1)*(c(2)-d(2))+c(1)*d(2))*(a(3)-d(3))+(a(1)-d(1)) &
        *(-c(3)*d(2)+b(3)*(-c(2)+d(2))+b(2)*(c(3)-d(3))+c(2)*d(3))+(a(2)-d(2))*(b(3)*(c(1)-d(1)) &
        +c(3)*d(1)-c(1)*d(3)+b(1)*(-c(3)+d(3))))/6
      totalvol1=totalvol1+tetvol(i)
    ENDDO
    WRITE(*,'(A,ES24.16)')'total volume:      ',totalvol1
    WRITE(*,'(A,ES24.16)')'equivalent radius: ',(3.0/4.0/pi*totalvol1)**(1.0/3.0)
  ENDSUBROUTINE computevol

END PROGRAM thrm2hammer