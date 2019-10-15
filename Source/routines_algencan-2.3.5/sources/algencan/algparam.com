C     COMMON SCALARS
      character * 6 hptype
      character * 2 innslvr
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,skipacc,
     +        sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,skipacc,sclsys,innercall,useustp,
     +                  innslvr,hptype
      save   /algparam/
