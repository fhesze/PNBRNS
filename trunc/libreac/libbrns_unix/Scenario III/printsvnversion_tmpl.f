       subroutine printSvnVersion()
           write(*,*) 'SVN Revision: $WCREV$'
           write(*,*) 'Committed on: $WCDATE$'
           write(*,*) 'Build time: $WCNOW$'
           write(*,*)
       end
