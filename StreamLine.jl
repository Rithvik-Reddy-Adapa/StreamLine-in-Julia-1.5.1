
using Interpolations

"""

Function Documentation

    #StreamLine(;x,y,z=nothing,u,v,w=nothing,xi,yi,zi=nothing,n,length,r_mag=false,r_u=false,r_v=false,r_w=false)

    #   Made on 10 / 10 / 2020 (dd / mm / yyyy)
    #   modules to use "using Interpolations"
    #    x,y,z,u,v,w = given vector data (u,v,w) for corresponding (x,y,z)
    #       x,y,z,u,v,w should be 3d array ( for 3d streamline )
    #       x,y,z should be grided as "x,y,z=[[i for i in x,j in y,k in z],[j for i in x,j in y,k in z],[k for i in x,j in y,k in z]]"
    #                                      OR
    #       x,y,u,v should be 2d array ( for 2d streamline )
    #       x,y should be grided as "x,y=[[i for i in x,j in y,k in z],[j for i in x,j in y,k in z]]"
    #   xi = initial x values as a row i.e. 1d row array
    #   yi = initial y values as a row i.e. 1d row array
    #   zi = initial z values as a row i.e. 1d row array
    #   n = total number of lines = total number of points (excluding initial point)
    #   length = length of each line = distance between consecutive points
    #   r_mag(return magnitude) = whether to return magnitude of vectors or not       |   all these
    #   r_u(return u) = whether to return i_cap(x component) of vectors or not        |   are useful
    #   r_v(return v) = whether to return j_cap(y component) of vectors or not        |   for
    #   r_w(return w) = whether to return k_cap(z component) of vectors or not        |   gradient plotting (i.e. useful to color the lines based on these values)
    #   returns in the following format:
    #       [xo,yo,zo,(mag),(i_cap),(j_cap),(k_cap)]
    #           mag is returned if r_mag==true
    #           i_cap is returned if r_u==true
    #           j_cap is returned if r_v==true
    #           k_cap is returned if r_w==true
    #       xo, yo,zo,mag,i_cap,j_cap,k_cap have same dimensions i.e. array shape = (size(xi),n+1)

"""
function StreamLine(;x,y,z=nothing,u,v,w=nothing,xi,yi,zi=nothing,n,length,r_mag=false,r_u=false,r_v=false,r_w=false)

        if zi!=nothing && size(z[1,1,:])[1]>1     # 3d StreamLine

                xo=zeros(size(xi)[1],n+1)
                yo=zeros(size(yi)[1],n+1)
                zo=zeros(size(zi)[1],n+1)

                if r_mag==true mag=zeros(size(xi)[1],n+1) else mag=nothing end
                if r_u==true i_cap=zeros(size(xi)[1],n+1) else i_cap=nothing end
                if r_v==true j_cap=zeros(size(xi)[1],n+1) else j_cap=nothing end
                if r_w==true k_cap=zeros(size(xi)[1],n+1) else k_cap=nothing end

                xo[:,1]=xi
                yo[:,1]=yi
                zo[:,1]=zi

                #=
                itp_a=interpolate((x[:,1,1],y[1,:,1],z[1,1,:]),u,Gridded(Linear()))
                itp_b=interpolate((x[:,1,1],y[1,:,1],z[1,1,:]),v,Gridded(Linear()))
                itp_c=interpolate((x[:,1,1],y[1,:,1],z[1,1,:]),w,Gridded(Linear()))

                etp_a=extrapolate(itp_a,Flat())
                etp_b=extrapolate(itp_b,Flat())
                etp_c=extrapolate(itp_c,Flat())
                =#

                #=
                etp_a=LinearInterpolation((x[:,1,1],y[1,:,1],z[1,1,:]),u,extrapolation_bc = Line())
                etp_b=LinearInterpolation((x[:,1,1],y[1,:,1],z[1,1,:]),v,extrapolation_bc = Line())
                etp_c=LinearInterpolation((x[:,1,1],y[1,:,1],z[1,1,:]),w,extrapolation_bc = Line())
                =#


                etp_a=LinearInterpolation((x[:,1,1],y[1,:,1],z[1,1,:]),u,extrapolation_bc = NaN)
                etp_b=LinearInterpolation((x[:,1,1],y[1,:,1],z[1,1,:]),v,extrapolation_bc = NaN)
                etp_c=LinearInterpolation((x[:,1,1],y[1,:,1],z[1,1,:]),w,extrapolation_bc = NaN)


                for j=1:size(xi)[1]

                        for i=1:n

                                a=etp_a(xo[j,i],yo[j,i],zo[j,i])
                                b=etp_b(xo[j,i],yo[j,i],zo[j,i])
                                c=etp_c(xo[j,i],yo[j,i],zo[j,i])

                                if r_mag==true mag[j,i]=(a^2+b^2+c^2)^(1/2) end
                                if r_u==true i_cap[j,i]=a end
                                if r_v==true j_cap[j,i]=b end
                                if r_w==true k_cap[j,i]=c end

                                xo[j,i+1]=xo[j,i]+(a*length)/((a^2+b^2+c^2)^(1/2))
                                yo[j,i+1]=yo[j,i]+(b*length)/((a^2+b^2+c^2)^(1/2))
                                zo[j,i+1]=zo[j,i]+(c*length)/((a^2+b^2+c^2)^(1/2))

                        end

                        a=etp_a(xo[j,n],yo[j,n],zo[j,n])
                        b=etp_b(xo[j,n],yo[j,n],zo[j,n])
                        c=etp_c(xo[j,n],yo[j,n],zo[j,n])

                        if r_mag==true mag[j,n]=(a^2+b^2+c^2)^(1/2) end
                        if r_u==true i_cap[j,n]=a end
                        if r_v==true j_cap[j,n]=b end
                        if r_w==true k_cap[j,n]=c end

                end

                ret=Any[xo,yo,zo,mag,i_cap,j_cap,k_cap]

                return ret

        else    # 2d StreamLine

                if z!=nothing && size(z[1,1,:])[1]==1 x=x[:,:,1]; y=y[:,:,1]; u=u[:,:,1]; v=v[:,:,1]; end

                xo=zeros(size(xi)[1],n+1)
                yo=zeros(size(yi)[1],n+1)

                if r_mag==true mag=zeros(size(xi)[1],n+1) else mag=nothing end
                if r_u==true i_cap=zeros(size(xi)[1],n+1) else i_cap=nothing end
                if r_v==true j_cap=zeros(size(xi)[1],n+1) else j_cap=nothing end

                xo[:,1]=xi
                yo[:,1]=yi

                #=
                itp_a=interpolate((x[:,1,1],y[1,:,1]),u,Gridded(Linear()))
                itp_b=interpolate((x[:,1,1],y[1,:,1]),v,Gridded(Linear()))

                etp_a=extrapolate(itp_a,Flat())
                etp_b=extrapolate(itp_b,Flat())
                =#

                #=
                etp_a=LinearInterpolation((x[:,1,1],y[1,:,1]),u,extrapolation_bc = Line())
                etp_b=LinearInterpolation((x[:,1,1],y[1,:,1]),v,extrapolation_bc = Line())
                =#


                etp_a=LinearInterpolation((x[:,1,1],y[1,:,1]),u,extrapolation_bc = NaN)
                etp_b=LinearInterpolation((x[:,1,1],y[1,:,1]),v,extrapolation_bc = NaN)


                for j=1:size(xi)[1]

                        for i=1:n

                                a=etp_a(xo[j,i],yo[j,i])
                                b=etp_b(xo[j,i],yo[j,i])

                                if r_mag==true mag[j,i]=(a^2+b^2)^(1/2) end
                                if r_u==true i_cap[j,i]=a end
                                if r_v==true j_cap[j,i]=b end

                                xo[j,i+1]=xo[j,i]+(a*length)/((a^2+b^2)^(1/2))
                                yo[j,i+1]=yo[j,i]+(b*length)/((a^2+b^2)^(1/2))

                        end

                        a=etp_a(xo[j,n],yo[j,n])
                        b=etp_b(xo[j,n],yo[j,n])

                        if r_mag==true mag[j,n]=(a^2+b^2)^(1/2) end
                        if r_u==true i_cap[j,n]=a end
                        if r_v==true j_cap[j,n]=b end

                end

                ret=Any[xo,yo,zo,mag,i_cap,j_cap]

                return ret

        end

end
