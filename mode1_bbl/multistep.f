C----------------------------------------------------------
C-Contains all functions that calculate BDF3 and AB3
C-coefficients for variable timestep.
C----------------------------------------------------------
C-AB3 coefficients
c
        function abvar1(dt1,dt2,dt3)

        dt123 = 2.*dt1 + 6.*dt2 + 3.*dt3
        dt23 = dt2 + dt3
        abvar1 = 1. + dt1*dt123/(6.*dt2*dt23)

        return
        end

c------------

        function abvar2(dt1,dt2,dt3)

        dt123 = 2.*dt1 + 3.*dt2 + 3.*dt3
        abvar2 = -dt1*dt123/(6.*dt2*dt3) 
    
        return
        end

c------------

        function abvar3(dt1,dt2,dt3)

        dt12 = 2.*dt1 + 3.*dt2
        dt23 = dt2 + dt3
        abvar3 = dt1*dt12/(6.*dt3*dt23)

        return
        end 
      
C-BDF3 coefficients
c
        function bdfvarv(dt1,dt2,dt3)

        dt12 = dt1+dt2
        dt123 = dt1+dt2+dt3 
        bdfvarv = 1. + dt1/dt12 + dt1/dt123

        return
        end

c-------------

        function bdfvar1(dt1,dt2,dt3)

        dt12 = dt1+dt2
        dt23 = dt2+dt3
        dt123 = dt1+dt2+dt3
        bdfvar1 = -1. - dt1/dt2 - dt1/dt123  
     >            -dt1*dt1/(dt123*dt2) 
     >            -dt1*dt1*dt12/(dt123*dt23*dt2) 
        bdfvar1 = -bdfvar1

        return
        end

c-------------

        function bdfvar2(dt1,dt2,dt3)

        dt12 = dt1+dt2
        dt23 = dt2+dt3
        dt123 = dt1+dt2+dt3
        bdfvar2 = dt1*dt1/(dt12*dt2) + dt1*dt1/(dt2*dt123)
     >           +dt1*dt1*dt12/(dt2*dt23*dt123)
     >           +dt1*dt1*dt12/(dt3*dt23*dt123)
        bdfvar2 = -bdfvar2 

        return
        end 

c--------------

        function bdfvar3(dt1,dt2,dt3)

        dt12 = dt1+dt2
        dt23 = dt2+dt3
        dt123 = dt1+dt2+dt3
        bdfvar3 = -dt1*dt1*dt12/(dt3*dt23*dt123)            
        bdfvar3 = -bdfvar3

        return
        end        
