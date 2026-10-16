program main
        implicit none
        character(len=4) :: info(8)
        info="O" 
        open(1000,file="info_connection")
        read(1000,*)
        read(1000,*)info(1:3)
        write(*,"(3a4)")info
        close(1000)
        


end




