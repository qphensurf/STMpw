! Normalizes the name of the elements
!      
      function norm_name( string ) 
      character(len=2)           :: string, norm_name
      integer                    :: i
      integer                    :: j
      norm_name = string
      j = iachar(string(1:1))
      if ( j >= iachar('a') .and. j <= iachar('z') ) then
          j = j + iachar('A') - iachar('a')
          norm_name(1:1) = achar(j)
      endif
      j = iachar(string(2:2))
      if ( j >= iachar('A') .and. j <= iachar('Z') ) then
          j = j + iachar('a') - iachar('A')
          norm_name(2:2) = achar(j)
      endif
      end function norm_name
