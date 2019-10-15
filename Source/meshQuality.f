c***********************************************************
      subroutine meshQuality2Dtris(nElTopo,topo,nVert,vertexes,crit)
c     compute the quality criterion for each element of topo
c     elements must be triangles (quadrilaterals not allowed yet)
c      USE IOUNIT
c      USE MAINMEM
      implicit none
      integer nElTopo! number of elements in topo
      integer nVert! number of vertexes in vertexes
      integer topo(nElTopo,3)! same format as lnods that stores the topology for the full mesh
      double precision vertexes(2,nVert)! coordinates of the vertexes in the mesh
      double precision crit(nElTopo)! value of the quality criterion for each element of topo
      integer el
      double precision vol,crtmp
      double precision M(2,2),Minv(2,2),W(2,2),Winv(2,2),J(2,2)
c----------------------------------------------------------
c     The quality criterion used is based on the squared condition number (extension of the work by Knupp)
      W(1,1)=1.0d0
      W(2,1)=0.0d0
      W(1,2)=0.5d0
      W(2,2)=sqrt(3.0d0)*0.5d0
      Winv(1,1)=1.0d0
      Winv(2,1)=0.0d0
      Winv(1,2)=-1.0d0/sqrt(3.0d0)
      Winv(2,2)=2.0d0/sqrt(3.0d0)
      do el=1,nElTopo
c        print*,"topo=",topo(el,:)
c        print*,"vert1=",vertexes(:,topo(el,1))
c        print*,"vert2=",vertexes(:,topo(el,2))
c        print*,"vert3=",vertexes(:,topo(el,3))
        
        J(1,1)=vertexes(1,topo(el,2))-vertexes(1,topo(el,1))
        J(2,1)=vertexes(2,topo(el,2))-vertexes(2,topo(el,1))
        J(1,2)=vertexes(1,topo(el,3))-vertexes(1,topo(el,1))
        J(2,2)=vertexes(2,topo(el,3))-vertexes(2,topo(el,1))
        M(1,1)=J(1,1)*Winv(1,1)+J(1,2)*Winv(2,1)
        M(1,2)=J(1,1)*Winv(1,2)+J(1,2)*Winv(2,2)
        M(2,1)=J(2,1)*Winv(1,1)+J(2,2)*Winv(2,1)
        M(2,2)=J(2,1)*Winv(1,2)+J(2,2)*Winv(2,2)
c        print*,"M=",M
        vol=M(1,1)*M(2,2)-M(1,2)*M(2,1)
        if (vol.gt.1.0d-8) then
          Minv(1,1)=M(2,2)/vol
          Minv(1,2)=-M(1,2)/vol
          Minv(2,1)=-M(2,1)/vol
          Minv(2,2)=M(1,1)/vol
c          print*,"Minv=",Minv
          M(1,1)=M(1,1)*M(1,1)+M(1,2)*M(1,2)+M(2,1)*M(2,1)+M(2,2)*M(2,2)! M(1,1) is re-used
          Minv(1,1)=Minv(1,1)*Minv(1,1)+Minv(1,2)*Minv(1,2)+
     &              Minv(2,1)*Minv(2,1)+Minv(2,2)*Minv(2,2)! Minv(1,1) is re-used
          crit(el)=0.25d0*M(1,1)*Minv(1,1)-1.0d0
c          print*,"cr=",crit(el)
        else
c          print*,"negative jacobian"
          crit(el)=1000.0d0! use a very large value meaning in inverted element
        endif
        if (crit(el).lt.0.0d0) crit(el)=0.0d0! for roundoffs
      enddo
      return
      end


c***********************************************************
      subroutine meshQuality2DtrisJF(nElTopo,topo,nVert,vertexes,crit)
c     compute the quality criterion for each element of topo
c     elements must be triangles (quadrilaterals not allowed yet)
c      USE IOUNIT
c      USE MAINMEM
      implicit none
      integer nElTopo! number of elements in topo
      integer nVert! number of vertexes in vertexes
      integer topo(nElTopo,3)! same format as lnods that stores the topology for the full mesh
      double precision vertexes(2,nVert)! coordinates of the vertexes in the mesh
      double precision crit(nElTopo)! value of the quality criterion for each element of topo
      integer el
      double precision vol,crtmp
      double precision M(2,2),Minv(2,2)
c----------------------------------------------------------
c     The quality criterion used is based on the squared condition number as defined by Knupp
c      note: I compute the SCN at every vertex and so aim at rectangular triangles, 
c            Knupp uses a projection to get equilateral triangles
      do el=1,nElTopo
c        print*,"topo=",topo(el,:)
c        print*,"vert1=",vertexes(:,topo(el,1))
c        print*,"vert2=",vertexes(:,topo(el,2))
c        print*,"vert3=",vertexes(:,topo(el,3))
        
        M(1,1)=vertexes(1,topo(el,2))-vertexes(1,topo(el,1))
        M(2,1)=vertexes(2,topo(el,2))-vertexes(2,topo(el,1))
        M(1,2)=vertexes(1,topo(el,3))-vertexes(1,topo(el,1))
        M(2,2)=vertexes(2,topo(el,3))-vertexes(2,topo(el,1))
c        print*,"M=",M
        vol=M(1,1)*M(2,2)-M(1,2)*M(2,1)
        if (vol.gt.1.0d-8) then
          Minv(1,1)=M(2,2)/vol
          Minv(1,2)=-M(1,2)/vol
          Minv(2,1)=-M(2,1)/vol
          Minv(2,2)=M(1,1)/vol
c          print*,"Minv=",Minv
          M(1,1)=M(1,1)*M(1,1)+M(1,2)*M(1,2)+M(2,1)*M(2,1)+M(2,2)*M(2,2)! M(1,1) is re-used
          Minv(1,1)=Minv(1,1)*Minv(1,1)+Minv(1,2)*Minv(1,2)+
     &              Minv(2,1)*Minv(2,1)+Minv(2,2)*Minv(2,2)! Minv(1,1) is re-used
          crit(el)=0.25d0*M(1,1)*Minv(1,1)-1.0d0
c          print*,"cr=",crit(el)
        else
c          print*,"negative jacobian"
          crit(el)=1000.0d0! use a very large value meaning in inverted element
        endif
        M(1,1)=vertexes(1,topo(el,3))-vertexes(1,topo(el,2))
        M(2,1)=vertexes(2,topo(el,3))-vertexes(2,topo(el,2))
        M(1,2)=vertexes(1,topo(el,1))-vertexes(1,topo(el,2))
        M(2,2)=vertexes(2,topo(el,1))-vertexes(2,topo(el,2))
        vol=M(1,1)*M(2,2)-M(1,2)*M(2,1)
        if (vol.gt.1.0d-8) then
          Minv(1,1)=M(2,2)/vol
          Minv(1,2)=-M(1,2)/vol
          Minv(2,1)=-M(2,1)/vol
          Minv(2,2)=M(1,1)/vol
          M(1,1)=M(1,1)*M(1,1)+M(1,2)*M(1,2)+M(2,1)*M(2,1)+M(2,2)*M(2,2)! M(1,1) is re-used
          Minv(1,1)=Minv(1,1)*Minv(1,1)+Minv(1,2)*Minv(1,2)+
     &              Minv(2,1)*Minv(2,1)+Minv(2,2)*Minv(2,2)! Minv(1,1) is re-used
          crtmp=0.25d0*M(1,1)*Minv(1,1)-1.0d0
        else
          crtmp=1000.0d0! use a very large value meaning in inverted element
        endif
        if (crtmp.lt.crit(el)) crit(el)=crtmp
        M(1,1)=vertexes(1,topo(el,1))-vertexes(1,topo(el,3))
        M(2,1)=vertexes(2,topo(el,1))-vertexes(2,topo(el,3))
        M(1,2)=vertexes(1,topo(el,2))-vertexes(1,topo(el,3))
        M(2,2)=vertexes(2,topo(el,2))-vertexes(2,topo(el,3))
        vol=M(1,1)*M(2,2)-M(1,2)*M(2,1)
        if (vol.gt.1.0d-8) then
          Minv(1,1)=M(2,2)/vol
          Minv(1,2)=-M(1,2)/vol
          Minv(2,1)=-M(2,1)/vol
          Minv(2,2)=M(1,1)/vol
          M(1,1)=M(1,1)*M(1,1)+M(1,2)*M(1,2)+M(2,1)*M(2,1)+M(2,2)*M(2,2)! M(1,1) is re-used
          Minv(1,1)=Minv(1,1)*Minv(1,1)+Minv(1,2)*Minv(1,2)+
     &              Minv(2,1)*Minv(2,1)+Minv(2,2)*Minv(2,2)! Minv(1,1) is re-used
          crtmp=0.25d0*M(1,1)*Minv(1,1)-1.0d0
        else
          crtmp=1000.0d0! use a very large value meaning in inverted element
        endif
        if (crtmp.lt.crit(el)) crit(el)=crtmp
        if (crit(el).lt.0.0d0) crit(el)=0.0d0! for roundoffs
      enddo
      return
      end


      subroutine meshQuality3Dtets(nElTopo,topo,nVert,vertexes,crit)
c     compute the quality criterion for each element of topo
c     elements must be tetras (hexas and pentas not allowed yet)
c      USE IOUNIT
c      USE MAINMEM
      implicit none
      integer nElTopo! number of elements in topo
      integer nVert! number of vertexes in vertexes
      integer topo(nElTopo,4)! same format as lnods that stores the topology for the full mesh
      double precision vertexes(nVert,3)! coordinates of the vertexes in the mesh
      double precision crit(nElTopo)! value of the quality criterion for each element of topo
      integer el
      double precision vol
      double precision M(3,3),Minv(3,3),W(3,3),Winv(3,3),J(3,3)
c----------------------------------------------------------
c     The quality criterion used is based on the squared condition number (extension of the work by Knupp)
      W(1,1)=1.0d0
      W(2,1)=0.0d0
      W(3,1)=0.0d0
      W(1,2)=0.5d0
      W(2,2)=sqrt(3.0d0)*0.5d0
      W(3,2)=0.0d0
      W(1,3)=0.5d0
      W(2,3)=sqrt(3.0d0)/6.0d0
      W(3,3)=sqrt(2.0d0/3.0d0)
      Winv(1,1)=1.0d0
      Winv(2,1)=0.0d0
      Winv(3,1)=0.0d0
      Winv(1,2)=-1.0d0/sqrt(3.0d0)
      Winv(2,2)=2.0d0/sqrt(3.0d0)
      Winv(3,2)=0.0d0
      Winv(1,3)=-1.0d0/sqrt(6.0d0)
      Winv(2,3)=-1.0d0/sqrt(6.0d0)
      Winv(3,3)=sqrt(3.0d0/2.0d0)
      do el=1,nElTopo
        J(1,1)=vertexes(topo(el,2),1)-vertexes(topo(el,1),1)
        J(2,1)=vertexes(topo(el,2),2)-vertexes(topo(el,1),2)
        J(3,1)=vertexes(topo(el,2),3)-vertexes(topo(el,1),3)
        J(1,2)=vertexes(topo(el,3),1)-vertexes(topo(el,1),1)
        J(2,2)=vertexes(topo(el,3),2)-vertexes(topo(el,1),2)
        J(3,2)=vertexes(topo(el,3),3)-vertexes(topo(el,1),3)
        J(1,3)=vertexes(topo(el,4),1)-vertexes(topo(el,1),1)
        J(2,3)=vertexes(topo(el,4),2)-vertexes(topo(el,1),2)
        J(3,3)=vertexes(topo(el,4),3)-vertexes(topo(el,1),3)
        M(1,1)=J(1,1)*Winv(1,1)+J(1,2)*Winv(2,1)+J(1,3)*Winv(3,1)
        M(2,1)=J(2,1)*Winv(1,1)+J(2,2)*Winv(2,1)+J(2,3)*Winv(3,1)
        M(3,1)=J(3,1)*Winv(1,1)+J(3,2)*Winv(2,1)+J(3,3)*Winv(3,1)
        M(1,2)=J(1,1)*Winv(1,2)+J(1,2)*Winv(2,2)+J(1,3)*Winv(3,2)
        M(2,2)=J(2,1)*Winv(1,2)+J(2,2)*Winv(2,2)+J(2,3)*Winv(3,2)
        M(3,2)=J(3,1)*Winv(1,2)+J(3,2)*Winv(2,2)+J(3,3)*Winv(3,2)
        M(1,3)=J(1,1)*Winv(1,3)+J(1,2)*Winv(2,3)+J(1,3)*Winv(3,3)
        M(2,3)=J(2,1)*Winv(1,3)+J(2,2)*Winv(2,3)+J(2,3)*Winv(3,3)
        M(3,3)=J(3,1)*Winv(1,3)+J(3,2)*Winv(2,3)+J(3,3)*Winv(3,3)
        vol=M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))-
     &      M(2,1)*(M(1,2)*M(3,3)-M(3,2)*M(1,3))+
     &      M(3,1)*(M(1,2)*M(2,3)-M(2,2)*M(1,3))
        if (vol.gt.1.0d-8) then
          Minv(1,1)=(M(2,2)*M(3,3)-M(2,3)*M(3,2))/vol
          Minv(1,2)=-(M(1,2)*M(3,3)-M(3,2)*M(1,3))/vol
          Minv(1,3)=(M(1,2)*M(2,3)-M(2,2)*M(1,3))/vol
          Minv(2,1)=-(M(2,1)*M(3,3)-M(3,1)*M(2,3))/vol
          Minv(2,2)=(M(1,1)*M(3,3)-M(3,1)*M(1,3))/vol
          Minv(2,3)=-(M(1,1)*M(2,3)-M(2,1)*M(1,3))/vol
          Minv(3,1)=(M(2,1)*M(3,2)-M(3,1)*M(2,2))/vol
          Minv(3,2)=-(M(1,1)*M(3,2)-M(3,1)*M(1,2))/vol
          Minv(3,3)=(M(1,1)*M(2,2)-M(1,2)*M(2,1))/vol
          M(1,1)=M(1,1)*M(1,1)+M(1,2)*M(1,2)+M(1,3)*M(1,3)+
     &           M(2,1)*M(2,1)+M(2,2)*M(2,2)+M(2,3)*M(2,3)+
     %           M(3,1)*M(3,1)+M(3,2)*M(3,2)+M(3,3)*M(3,3)! M(1,1) is re-used
          Minv(1,1)=Minv(1,1)*Minv(1,1)+Minv(1,2)*Minv(1,2)+
     &              Minv(1,3)*Minv(1,3)+Minv(2,1)*Minv(2,1)+
     &              Minv(2,2)*Minv(2,2)+Minv(2,3)*Minv(2,3)+
     &              Minv(3,1)*Minv(3,1)+Minv(3,2)*Minv(3,2)+
     &              Minv(3,3)*Minv(3,3)! Minv(1,1) is re-used
          crit(el)=0.1111111111111d0*M(1,1)*Minv(1,1)-1.0d0
        else
          crit(el)=1000.0d0! use a very large value meaning in inverted element
        endif
        if (crit(el).lt.0.0d0) crit(el)=0.0d0! for roundoffs
      enddo
      return
      end

