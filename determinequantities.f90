    module determinequantities
    
      use iso_fortran_env, only: int64
! Code for reading a NORMWAVECAR function
!  coded by A. Arnau and N. Lorente
!
!last edit
!    05/04/20
     
contains

      subroutine determine_quantities( Name, Number_of_SPIN, Total_Bands, &
     &           Number_of_KPOINTS)
      implicit none

      INTEGER, PARAMETER :: q =SELECTED_REAL_KIND(10)
      INTEGER, PARAMETER :: qs=SELECTED_REAL_KIND(5)
        INTEGER, PARAMETER :: ICMPLX=16
        integer(int64) :: wavecar_size, expected_size
        integer(int64) :: recl64, nrecords64
      INTEGER, PARAMETER ::  IU0=6

      REAL(q)  RDUM, RISPIN, norma, ENMAXF
      REAL(q)  RTAG, RKPTSF, RBANDF, RNPL
      REAL(q)  A(3,3)
      character( len = 60) Name
      INTEGER IDUM, I, J, K, ISP, ISPIN, NKPTS, NBANDS
      integer IREC,Total_Bands, Number_of_SPIN, Number_of_KPOINTS
! first reopen with assumed (wrong) record length ICMPLX
      OPEN(UNIT=12,FILE=Name,ACCESS='DIRECT', &
                         FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=ICMPLX)

! the first record contains the record length, get it ...
      RDUM=0._q
      READ(12,REC=1,ERR=17421) RDUM,RISPIN ; IDUM=NINT(RDUM)
        IF (IDUM<=0) IDUM=ICMPLX  ! error reading RECL
      GOTO 17422
17421 CONTINUE
      IDUM=ICMPLX
17422 CONTINUE
      close (12)

! reopen with correct record length (clumsy all that, I know ...)
      OPEN(UNIT=12,FILE=Name,ACCESS='DIRECT', &
                      FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IDUM)



! header dimensions and alloctations

        RTAG=0
        READ(12,REC=1,ERR=200) RDUM,RISPIN,RTAG

        if (RTAG/=45200) then
                write(IU0,*) 'Please, change complex(qs) to complex(q)'
                write(IU0,*) 'and recompile readwf.f90'
                write(IU0,*) 'Thank you for your cooperation.'
                write(IU0,*) '(do not forget to remove this if)'
                stop
        endif

        ISPIN=NINT(RISPIN)
        Number_of_SPIN = ISPIN

        IREC=2
        READ(12,REC=2,ERR=200) RKPTSF,RBANDF,ENMAXF, &
             &                         ((A(I,J),I=1,3),J=1,3)

        IREC=2
        NKPTS=NINT(RKPTSF)
        NBANDS=NINT(RBANDF)
        Number_of_KPOINTS = NKPTS
        Total_Bands = NBANDS

! Validate WAVECAR size using explicitly 64-bit arithmetic.
          inquire(unit=12,size=wavecar_size)
          recl64=int(IDUM,int64)
          nrecords64=2_int64+int(ISPIN,int64)*int(NKPTS,int64)* &
       &             (int(NBANDS,int64)+1_int64)
          expected_size=recl64*nrecords64
          if(wavecar_size.ne.expected_size) then
            write(IU0,'(A,I0)') 'WAVECAR size       = ',wavecar_size
            write(IU0,'(A,I0)') 'Expected size      = ',expected_size
            write(IU0,'(A,I0)') 'WAVECAR RECL bytes = ',recl64
            write(IU0,'(A,I0)') 'Number of records  = ',nrecords64
            write(IU0,*) 'WAVECAR size does not match its header'
            stop
          endif
	close (12)

	return
200   write (IU0,*) 'WAVECAR is corrupt (check its size!)'
      stop
230   write (IU0,*) 'Error while reading Eigenvalues in WAVECAR'
      write (IU0,*) K,ISP
      stop
240   write (IU0,*) 'Error while reading EigenVECTEURS in WAVECAR'
      write (IU0,*) J,K,ISP
      stop

	end subroutine

    end module determinequantities



