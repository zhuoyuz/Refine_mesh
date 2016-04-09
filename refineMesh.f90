program refine_mesh
   implicit none
   integer :: i,j,k,ifile
   integer :: ie,iAdd,iElem
   integer :: iNode, jNode
   integer :: nNodes=1302,nElems=2464
   integer :: nNodesMax,nElemsMax
   integer :: nAdd
   integer :: nfile = 1001
   double precision, allocatable :: x(:),y(:),z(:)
   double precision, allocatable :: xFine(:),yFine(:),zFine(:)
   double precision :: x1, y1, z1
   double precision :: x2, y2, z2
   double precision :: x3, y3, z3
   double precision :: eps = 1.e-6
   integer, allocatable :: id(:) !Node index
   integer, allocatable :: tri(:,:)
   integer, allocatable :: triFine(:,:)
   integer, allocatable :: elem(:,:)
   logical, allocatable :: visited(:), discard(:)
   character*30 :: fname

   nNodesMax = nNodes + nElems*3
   nElemsMax = nElems*4
   nAdd = nElems*3

   allocate(x(nNodes),y(nNodes),z(nNodes))
   allocate(xFine(nNodesMax),yFine(nNodesMax),zFine(nNodesMax))

   allocate(id(nNodesMax))
   allocate(tri(3,nElems))
   allocate(triFine(3,nElemsMax))

   allocate(elem(3,nNodes+1:nNodesMax))
   allocate(visited(nNodes+1:nNodesMax), discard(nNodes+1:nNodesMax))

   do iNode = 1, nNodesMax
      id(iNode) = iNode
   enddo

   ! Read Coarse Data File
   ifile = 10000
   open(ifile,file='marker.dat',status='unknown')
   do i = 1, 4
      read(ifile,*)
   enddo
   do j = 1, nfile, 1
      do i = 1, 5
         read(ifile,*)
      enddo
      do i = 1, nNodes
         read(ifile,*) x(i),y(i),z(i)
      enddo
      do ie = 1, nElems
         read(ifile,*) tri(1,ie), tri(2,ie), tri(3,ie)
      enddo

      xFine(1:nNodes) = x(1:nNodes)
      yFine(1:nNodes) = y(1:nNodes)
      zFine(1:nNodes) = z(1:nNodes)

      iAdd = nNodes
      iElem = 0
      do ie =1, nElems
         x1 = (x(tri(1,ie))+x(tri(2,ie)))/2
         y1 = (y(tri(1,ie))+y(tri(2,ie)))/2
         z1 = (z(tri(1,ie))+z(tri(2,ie)))/2
         x2 = (x(tri(2,ie))+x(tri(3,ie)))/2
         y2 = (y(tri(2,ie))+y(tri(3,ie)))/2
         z2 = (z(tri(2,ie))+z(tri(3,ie)))/2
         x3 = (x(tri(3,ie))+x(tri(1,ie)))/2
         y3 = (y(tri(3,ie))+y(tri(1,ie)))/2
         z3 = (z(tri(3,ie))+z(tri(1,ie)))/2

         xFine(iAdd+1) = x1
         yFine(iAdd+1) = y1
         zFine(iAdd+1) = z1
         xFine(iAdd+2) = x2
         yFine(iAdd+2) = y2
         zFine(iAdd+2) = z2
         xFine(iAdd+3) = x3
         yFine(iAdd+3) = y3
         zFine(iAdd+3) = z3

         triFine(1,iElem+1) = tri(1,ie)
         triFine(2,iElem+1) = iAdd + 1
         triFine(3,iElem+1) = iAdd + 3

         triFine(1,iElem+2) = iAdd + 1
         triFine(2,iElem+2) = tri(2,ie)
         triFine(3,iElem+2) = iAdd + 2

         triFine(1,iElem+3) = iAdd + 2
         triFine(2,iElem+3) = tri(3,ie)
         triFine(3,iElem+3) = iAdd + 3

         triFine(1,iElem+4) = iAdd + 1
         triFine(2,iElem+4) = iAdd + 2
         triFine(3,iElem+4) = iAdd + 3

         ! Each newly added node has 3 neighboring elems
         elem(1,iAdd+1) = iElem + 1
         elem(2,iAdd+1) = iElem + 4
         elem(3,iAdd+1) = iElem + 2

         elem(1,iAdd+2) = iElem + 2
         elem(2,iAdd+2) = iElem + 4
         elem(3,iAdd+2) = iElem + 3

         elem(1,iAdd+3) = iElem + 3
         elem(2,iAdd+3) = iElem + 4
         elem(3,iAdd+3) = iElem + 1

         iAdd = iAdd + 3
         iElem = iElem + 4
      enddo

      print *, iAdd, iElem

      !do iNode = 1, nNodesMax
      !   write(101,'(3f15.10)') xFine(iNode), yFine(iNode), zFine(iNode)
      !enddo

      ! Combine nodes with same coordinates
      visited = .false.
      discard = .false.
      do iNode = nNodes+1, nNodesMax
         if (.not.visited(iNode)) then
            visited(iNode) = .true.
            do jNode = iNode+1, nNodesMax
               if (.not.visited(jNode)) then
                  ! Combine iNode and jNode
                  if (abs(xFine(iNode)-xFine(jNode))<eps .and. &
                     abs(yFine(iNode)-yFine(jNode))<eps .and. &
                     abs(zFine(iNode)-zFine(jNode))<eps) then
                     visited(jNode) = .true.
                     discard(jNode) = .true.
                     !write(103,'(3f15.10)') xFine(iNode),yFine(iNode),zFine(iNode)
                     !write(103,'(3f15.10)') xFine(jNode),yFine(jNode),zFine(jNode)
                     !write(103,*) iNode, jNode
                     !write(103,*) 
                     ! Element # should not change
                     if (triFine(1,elem(1,jNode)) == jNode) &
                        triFine(1,elem(1,jNode)) = iNode
                     if (triFine(2,elem(1,jNode)) == jNode) &
                        triFine(2,elem(1,jNode)) = iNode
                     if (triFine(3,elem(1,jNode)) == jNode) &
                        triFine(3,elem(1,jNode)) = iNode
                     if (triFine(1,elem(2,jNode)) == jNode) &
                        triFine(1,elem(2,jNode)) = iNode
                     if (triFine(2,elem(2,jNode)) == jNode) &
                        triFine(2,elem(2,jNode)) = iNode
                     if (triFine(3,elem(2,jNode)) == jNode) &
                        triFine(3,elem(2,jNode)) = iNode
                     if (triFine(1,elem(3,jNode)) == jNode) &
                        triFine(1,elem(3,jNode)) = iNode
                     if (triFine(2,elem(3,jNode)) == jNode) &
                        triFine(2,elem(3,jNode)) = iNode
                     if (triFine(3,elem(3,jNode)) == jNode) &
                        triFine(3,elem(3,jNode)) = iNode
                  endif
               endif
            enddo
         endif
      enddo

      !do iNode = 1, nAdd+nNodes
      !   if (.not.discard(iNode)) then
      !      write(104,'(3f15.10)') xFine(iNode), yFine(iNode), zFine(iNode)
      !   endif
      !enddo

      nAdd = nElems*3
      do iNode = nNodesMax, nNodes+1, -1
         if (discard(iNode)) then
            !print *, "dicarded nodes", iNode
            ! Shift all right nodes one step to the left
            do jNode = iNode+1, nAdd+nNodes
               id(jNode) = jNode - 1
               do ie = 1, nElemsMax
                  if (triFine(1,ie)==jNode) triFine(1,ie) = jNode-1
                  if (triFine(2,ie)==jNode) triFine(2,ie) = jNode-1
                  if (triFine(3,ie)==jNode) triFine(3,ie) = jNode-1
                  xFine(jNode-1) = xFine(jNode)
                  yFine(jNode-1) = yFine(jNode)
                  zFine(jNode-1) = zFine(jNode)
               enddo
            enddo
            nAdd = nAdd - 1
         endif
      enddo

      !do iNode = 1, nAdd+nNodes
      !   write(102,'(3f15.10)') xFine(iNode), yFine(iNode), zFine(iNode)
      !enddo

      print *, nAdd+nNodes

      ! output
      write(fname,'(A,I5.5,A)') 'fineMesh.',j,'.dat'
      open(ifile+j,file=fname,status='replace')
      write(ifile+j,*) 'Variables = "X","Y","Z"'
      write(ifile+j,*) 'ZONE N=', nNodes+nAdd, ',E=', nElemsMax, ', ZONETYPE=FETRIANGLE'
      write(ifile+j,*) 'DATAPACKING=POINT'
      do i = 1, nNodes+nAdd
         write(ifile+j,'(3f15.10)') xFine(i),yFine(i),zFine(i)
      enddo
      do i = 1, nElemsMax
         write(ifile+j,'(3i7)') triFine(1,i),triFine(2,i),triFine(3,i)
      enddo
      close(ifile+j)

   enddo

   close(ifile)

end
