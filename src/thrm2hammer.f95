PROGRAM thrm2hammer
  IMPLICIT NONE
  !globals are up here as is the driver
  INTEGER(8) :: num_points,num_tets
  CHARACTER(64) :: mesh_file
  INTEGER(8),ALLOCATABLE :: tetgeom(:,:)
  REAL(8),ALLOCATABLE :: point(:,:)
  REAL(8),PARAMETER :: PI=4.D0*DATAN(1.D0)
  INTEGER(8),PARAMETER::in_unit=20,out_unit=30

  CALL readin_cmdline()

  CALL readin_mesh()

  CALL computevol()

  CALL outputthrm()

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
    OPEN(UNIT = in_unit, FILE = mesh_file, ACTION = "read", STATUS='OLD', IOSTAT=tempint)
    IF(tempint .NE. 0)STOP 'error opening input mesh file'

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
    REAL(8) :: totalvol1,tetvol
    REAL(8) :: a(3),b(3),c(3),d(3)

    totalvol1=0
    !compute the tet volumes
    DO i=1,num_tets
      a(:)=point(tetgeom(i,1),:)
      b(:)=point(tetgeom(i,2),:)
      c(:)=point(tetgeom(i,3),:)
      d(:)=point(tetgeom(i,4),:)
      tetvol=ABS((-c(2)*d(1)+b(2)*(-c(1)+d(1))+b(1)*(c(2)-d(2))+c(1)*d(2))*(a(3)-d(3))+(a(1)-d(1)) &
        *(-c(3)*d(2)+b(3)*(-c(2)+d(2))+b(2)*(c(3)-d(3))+c(2)*d(3))+(a(2)-d(2))*(b(3)*(c(1)-d(1)) &
        +c(3)*d(1)-c(1)*d(3)+b(1)*(-c(3)+d(3))))/6
      totalvol1=totalvol1+tetvol
    ENDDO
    WRITE(*,'(A,ES24.16)')'total volume:      ',totalvol1
    WRITE(*,'(A,ES24.16)')'equivalent radius: ',(3.0/4.0/pi*totalvol1)**(1.0/3.0)
  ENDSUBROUTINE computevol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE outputthrm()
    CHARACTER(64) :: outfile
    INTEGER(8) :: i,tempint
    !open the output mesh file
    WRITE(outfile,'(2A)')TRIM(ADJUSTL(mesh_file)),'.xml'
    OPEN(UNIT = out_unit, FILE = outfile, ACTION = "WRITE", STATUS="REPLACE", IOSTAT=tempint)
    IF(tempint .NE. 0)STOP 'error opening output mesh file'

    !initial wrapping stuff
    WRITE(out_unit,'(A)')"<?xml version='1.0' encoding='UTF-8'?>"
    WRITE(out_unit,*)
    WRITE(out_unit,'(A)')"<hammer>"
    WRITE(out_unit,'(A)')"  <external_meshes>"
    WRITE(out_unit,'(A,I0,A,I0,A)')'    <unstructured_tetrahedral_mesh name="my_tet_mesh" bins="', &
      num_tets,'" k="',num_tets,'" adjacency_info="false">'
    WRITE(out_unit,'(A)')"      <vertices>"
    !write out the points
    DO i=1,num_points
      WRITE(out_unit,'(A,3ES24.16)')'       ',point(i,:)
    ENDDO
    WRITE(out_unit,'(A)')"      </vertices>"
    WRITE(out_unit,'(A)')"      <elements>"
    !write out each tet
    DO i=1,num_tets
      WRITE(out_unit,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)')'        <e id="',i-1,'">',tetgeom(i,1)-1,' ', &
        tetgeom(i,2)-1,' ',tetgeom(i,3)-1,' ',tetgeom(i,4)-1,'</e>'
    ENDDO
    !closing wrapups
    WRITE(out_unit,'(A)')"      </elements>"
    WRITE(out_unit,'(A)')'    </unstructured_tetrahedral_mesh>'
    WRITE(out_unit,'(A)')"  </external_meshes>"
    WRITE(out_unit,'(A)')"</hammer>"

    WRITE(*,'(A)')'Program completed without error'
    WRITE(*,'(3A)')'Output mesh written to "',TRIM(ADJUSTL(outfile)),'"'

    CLOSE(out_unit)
  ENDSUBROUTINE outputthrm

END PROGRAM thrm2hammer