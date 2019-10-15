C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

