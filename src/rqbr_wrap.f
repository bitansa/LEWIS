C     Wrapper of the routine rqbr of the file rqbr.f.
C
C     The C++ code passes the flag lci1 as an INTEGER. This wrapper
C     turns it into the LOGICAL argument that rqbr expects and then
C     calls rqbr with all other arguments unchanged. The C++ declaration
C     then uses only INTEGER and DOUBLE PRECISION arguments. These types
C     match the Fortran code exactly.
      subroutine lwrqbr(m,nn,m5,n3,n4,a,b,t,toler,ift,x,e,s,wa,wb,
     *nsol,ndsol,sol,dsol,lsol,h,qn,cutoff,ci,tnmat,big,ilci1)
      integer m,nn,m5,n3,n4,ift,nsol,ndsol,lsol,ilci1
      integer s(m),h(nn,nsol)
      double precision a(m,nn),b(m),t,toler,x(nn),e(m)
      double precision wa(m5,n4),wb(m),sol(n3,nsol),dsol(m,ndsol)
      double precision qn(nn),cutoff,ci(4,nn),tnmat(4,nn),big
      logical lci1
C     Convert the INTEGER flag to a LOGICAL value.
      lci1 = ilci1 .ne. 0
C     Call the Barrodale-Roberts solver.
      call rqbr(m,nn,m5,n3,n4,a,b,t,toler,ift,x,e,s,wa,wb,nsol,ndsol,
     *sol,dsol,lsol,h,qn,cutoff,ci,tnmat,big,lci1)
      return
      end
