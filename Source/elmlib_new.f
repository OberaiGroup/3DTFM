c***********************************************************
      subroutine elmlib (ielem,itask,pelem,nnode,est,efo,eln,xel,
     $     elemd,uel,uel_meas)
c     element library

c     itask = 1  :  initialize element parameters
c             2  :  read material properties
c             3  :  compute element tangent stiffness
c             4  :  compute element internal force vector
c             5  :  compute the elemental RHS/residual contribution from a external source
c             6  :  build the elemental RHS/residual for the dual problem
c             7  :  compute the elemental contributions to the objective function and its gradient
c     Tasks 8 to 17 are not used (2010-03-05)
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ielem, itask
      integer nnode! localized information
      integer ltype
      double precision pelem(*)! localized information
      double precision est(mevab,mevab)! elemental tangent stiffness matrix
      double precision efo(mevab)! elemental rhs/residual
      integer eln(mnode)! nodes number in the element
      double precision xel(ndime,*)! coordinates of the nodes in the element
      double precision elemd(nset_nodal,mnode)! values of the material parameter in the element
      double precision uel(mdofn,mnode)! value of the displacements at the nodes of the element
      double precision uel_meas(mdofn,mnode)! value of the measured disp. at the nodes of the element
c----------------------------------------------------------

      ltype = lrefn(ielem,1)

c     LIST OF ELEMENTS
      if (ltype.eq.66) then
        call elem66 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.37) then
        call elem37 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.77) then
        call elem77 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.606) then
        call elem606 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.611) then
        call elem611 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.618) then
        call elem618 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.607) then
        call elem607 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.608) then
        call elem608 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.505) then
        call elem505 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.707) then
        call elem707 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.306) then
        call elem306 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      elseif (ltype.eq.305) then
        call elem305 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      else if (ltype.eq.307) then
        call elem307 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      else if (ltype.eq.31) then
        call elem31 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      else if (ltype.eq.32) then
        call elem32 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      else if (ltype.eq.309) then
        call elem309 (ielem,itask,pelem,nnode,est,efo,eln,xel,elemd,uel,
     $        uel_meas)
      else
        write(iwrit,900) ltype
        call termin(6)
      endif
      return

900   format(//2x,'** fatal error ** element type number',
     +        ' = ',i4,' not in element library')
      end

