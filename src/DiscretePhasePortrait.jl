module DiscretePhasePortrait

using Plots, Contour, LinearAlgebra
    
export phaseportrait

"""
	phaseportrait(F, G, limits; <keyword arguments>)

Plot a phase portrait for a discrete 2D system.

A phase portrait of the discrete-time system

```math
	x_{t+1} =  F(x_t, y_t),     y_{t+1} = G(x_t, y_t)
```

is drawn showing isoclines at particular directions and their preimages for the rectangular 
region of the plane defined by `limits=[xmin xmax; ymin ymax]`.  Labels and a legend are used 
to indicate which regions map to others and the direction of motion in each region.

# Keyword Arguments

  - `fixpt`: Approximation `(x0,y0)` of a fixed point of interest.
  - `Jac`: Function of `x` and `y` returning the Jacobian matrix of the system at a point `(x,y)`.
  - `directions=[0,π/2]`: Directions of the isoclines to use if `fixpt` or `Jac` are not given.
    !!! info
        If both `fixpt` and `Jac` are specified, the actual fixed point is determined using
        Newton's method and `Jac` starting from the input `fixpt`, then the eigenvalues and
        eigenvectors of the Jacobian at the actual fixed point are calculated and the
        directions of the isoclines are set to the eigenvector directions.  The return
        value is `(actualfixedpoint,eigenvalues,eigenvectors,evdirections)`.  If either
        `fixpt` or `Jac` are *not* specified then the directions of the isoclines are
        determined by the values in `directions`, and the return value is `nothing`.

  - `n::Int=1`: Number of iterations of the function f=[F,G] applied.
  - `labelpoints`: User defined points at which to apply region labeling.  

    This may be an array or tuple of `(x,y)` points.  If `labelpoints` is not specified then 
    by default, four points on the main plot near the four corners are labelled.

  - `showdet0::Bool=true`: Plot the det(Jac)=0 curve and its image. `Jac` must also be given.
  - `detlimits=limits`: Limits of the form `[dxmin dxmax; dymin dymax]` to use for the image of the 
    det(Jac)=0 curve.
    The points within `detlimits` satisfying det(Jac)=0 are computed and the image of all such
    points are plotted if they lie within `limits`. 
    !!! tip 
        Often, to see more of the image curve it is 
        necessary to make `detlimits` a larger rectangle than `limits`.

  - `showrange::Bool=true`: Gray shade the image of all points within `limits`.
  - `resolution::Int=500`: Grid resolution.
  - `plotsize=(800,600)`: Set the size of the plot in pixels.
"""
function phaseportrait(F, G, limits; n=1, directions=[0,π/2], fixpt=nothing, Jac=nothing, 
		showdet0=true, showrange=true, detlimits=nothing, labelpoints=nothing, resolution=500, 
		plotsize=(800,600))
	
	# If the limits are specified as a 4-element vector [xmin,xmax,ymin,ymax], reorder to 
	# a 2x2 matrix with x limits in the first row
	if size(limits) != (2,2)
		limits = [limits[1] limits[2];limits[3] limits[4]]
	end
	# Set up ranges of x and y values.
	x = range(limits[1,1], limits[1,2], length=resolution)
	y = range(limits[2,1], limits[2,2], length=resolution)'
    #  Use detlimits for plotting the image of the det J=0 line if given
	if showdet0 && !isnothing(detlimits) && !isnothing(Jac)
		if size(detlimits) != (2,2)
			detlimits = [detlimits[1] detlimits[2]; detlimits[3] detlimits[4]]
		end
		x2 = range(detlimits[1,1], detlimits[1,2], length=resolution)
		y2 = range(detlimits[2,1], detlimits[2,2], length=resolution)'
	else
		x2 = nothing
		y2 = nothing
	end

    # Define iterated functions.
	FG(x,y) = fun_iter(F, G, n, x, y)  # FG will return the tuple (F^n(x,y), G^n(x,y))

	if !isnothing(Jac)
		FG_Jac(x,y) = Jac_iter(F, G, Jac, n, x, y)
	end

    # If a fixed point estimate and the Jacobian are provided then use Newton's Method 
	# to find the fixed point and use the eigenvectors as the directions, θ.
	if !isnothing(fixpt) && !isnothing(Jac)
		actualfixpt = copy(fixpt)
		newton!(actualfixpt, X->collect(FG(X[1],X[2])) .- X, X->FG_Jac(X[1],X[2]) - I, 
				maximum(limits[:,2] .- limits[:,1])*1.0e-8)
		(evalues,evectors,direc) = getdirections(FG_Jac(actualfixpt[1],actualfixpt[2]))
		returnvalues = (actualfixpt,evalues,evectors,direc)
	else
		direc = directions
		returnvalues = nothing
    end
	
	# Size of plot in pixels
	(X_max, Y_max) = plotsize

	# colour selections for the two directions and det=0
    col = [:blue,:red]
	detcolor = :green

    # Initialize the plot for isoclines, preimages and θ angled grid lines
	p_1 = plot(xlabel="x", ylabel="y", xlims=limits[1,:], ylims=limits[2,:], grid=false, 
			   size=plotsize, legend=false, clf=true, framestyle=:box, widen=true)

    #This is the region where all the points in the plane maps
	if showrange
		forwardmap!(p_1,FG,x,y)
	end

    # generate grid lines for each θ using the function grid
    grid_lines = grid!(p_1,direc, limits, X_max, Y_max, col)
    
    # Now we plot isoclines and preimages for each θ.
    curves = iso_pre!(p_1,FG, direc, x, y, col)

	#ploting isocines and preimages at the edges
	x_val= limits[1,:]
	y_val= limits[2,:]
	for (i,θ) in enumerate(direc)
		#Left edge
	    C_left = (FG(x_val[1], y)[1] .- x_val[1]) .* sin(θ) .- (FG(x_val[1], y)[2] .- y) .* cos(θ)
		C_left_pre = ((FG(FG(x_val[1], y)[1], FG(x_val[1], y)[2])[1] .- FG(x_val[1], y)[1]) .* sin(θ) .- 
		             (FG(FG(x_val[1], y)[1], FG(x_val[1], y)[2])[2] .- FG(x_val[1], y)[2]) .* cos(θ))
		left = iso_pre_edge!(p_1, y, C_left, true, x_val[1], :solid, col[i])
		left_pre = iso_pre_edge!(p_1, y, C_left_pre, true, x_val[1],:dash, col[i])
	    #Right edge
	    C_right = (FG(x_val[2], y)[1] .- x_val[2]) .* sin(θ) .- (FG(x_val[2], y)[2] .- y) .* cos(θ)
		C_right_pre = ((FG(FG(x_val[2], y)[1], FG(x_val[2], y)[2])[1] .- FG(x_val[2], y)[1]) .* sin(θ) .- 
		              (FG(FG(x_val[2], y)[1], FG(x_val[2], y)[2])[2] .- FG(x_val[2], y)[2]) .* cos(θ))
	    right = iso_pre_edge!(p_1, y, C_right, true, x_val[2],:solid, col[i])
		right_pre = iso_pre_edge!(p_1, y, C_right_pre, true, x_val[2],:dash, col[i])
	    #Bottom edge
	    C_bottom = (FG(x, y_val[1])[1] .- x) .* sin(θ) .- (FG(x, y_val[1])[2] .- y_val[1]) .* cos(θ)
		C_bottom_pre = ((FG(FG(x, y_val[1])[1], FG(x, y_val[1])[2])[1] .- FG(x, y_val[1])[1]) .* sin(θ) .- 
		               (FG(FG(x, y_val[1])[1], FG(x, y_val[1])[2])[2] .- FG(x, y_val[1])[2]) .* cos(θ))
	    bottom = iso_pre_edge!(p_1, x, C_bottom, false, y_val[1],:solid, col[i])
		bottom_pre = iso_pre_edge!(p_1, x, C_bottom_pre, false, y_val[1],:dash, col[i])
	    #Top edge
	    C_top = (FG(x, y_val[2])[1] .- x) .* sin(θ) .- (FG(x, y_val[2])[2] .- y_val[2]) .* cos(θ)
		C_top_pre = ((FG(FG(x, y_val[2])[1], FG(x, y_val[2])[2])[1] .- FG(x, y_val[2])[1]) .* sin(θ) .- 
		            (FG(FG(x, y_val[2])[1], FG(x, y_val[2])[2])[2] .- FG(x, y_val[2])[2]) .* cos(θ))
	    top = iso_pre_edge!(p_1, x, C_top, false, y_val[2],:solid, col[i])
		top_pre = iso_pre_edge!(p_1, x, C_top_pre, false, y_val[2],:dash, col[i])
	end

    # Plot where the determinant is zero and the image of where the determinant is zero.
	if showdet0 && !isnothing(Jac)
		detplotter!(p_1,x,y,FG,FG_Jac,detcolor,x2,y2)
    end

	# Annotate points in the four corners of the plot
	plot_annotate!(p_1,limits,FG,direc,labelpoints)

	# create the legend plot
	p_2 = legendplot(limits,X_max,Y_max,direc,FG,col)

    # Create a layout with the main plot and the legend plot, and plot them
    L = @layout [a{0.85w} [b{0.2h}; c]]
	p = plot(p_1, p_2, layout = L, size=(plotsize[1]/0.85,plotsize[2]))
    gui(p)
	return returnvalues
end

# Compute the function [F,G] iterated n times, works for vector x, y values
# The return is a tuple of F^n values and G^n values.
function fun_iter(F,G,n,x,y)
    val_F = F.(x,y)
    val_G = G.(x,y)
    for i=2:n
        (val_F, val_G) = (F.(val_F, val_G), G.(val_F, val_G))
    end
    return (val_F, val_G)
end

# Compute the Jacobian of the n iterated mapping, only needed for scalar x,y values.
# Returns a 2x2 matrix.
function Jac_iter(F, G, Jac, n, x, y)
    J = Jac(x,y)
    val_F = x
    val_G = y
    for i = 2:n 
        (val_F, val_G) = (F(val_F, val_G), G(val_F, val_G))
        J = Jac(val_F, val_G) * J
    end
    return J 
end

function newton!(x, f, Df, Tol)
    # Compute a zero of f using Newton's method.  Df is the Jacobian, and Tol is the
	# tolerance.  On input x is the initial guess.  On output x is the zero.
	maxiter = 15
	s = f(x)  # If the initial guess is good enough, no iterations will be performed.
    k = 0
	while norm(s) > Tol && (k += 1) <= maxiter
		fx, Jx = f(x), Df(x)
        s = Jx\fx
        x .-= s
    end
    if k > maxiter
        error("Newton's method did not converge within $maxiter iterations.  "*
			  "Try a better initial guess for the fixed point.")
    end
end

function getdirections(M)
    # Get the directions of the eigenvectors of M.
	eigvals, eigvecs = eigen(M)  
    if all(isreal, eigvals)
        # If eigenvalues are real, use x and y components of eigenvectors to compute directions.
        direc = [atan(eigvecs[2,1]/eigvecs[1,1]), atan(eigvecs[2,2]/eigvecs[1,2])]
    else
        # If eigenvalues are complex, use real and imaginary parts of eigenvectors to compute directions.
        # Since one eigen vector is the complex conjugate of the other, no need to consider both.
		direc = [atan(imag(eigvecs[2,1])/imag(eigvecs[1,1])), atan(real(eigvecs[2,1])/real(eigvecs[1,1]))]
    end 
    for i in 1:2
        direc[i]=mod(direc[i],pi)
    end
	if direc[1]>direc[2]
		direc[1],direc[2]=direc[2],direc[1]
	end
	return (eigvals,eigvecs,direc)
end

#  Data coordinates are in the rectangle 
#        limits[1,1] <= x <= limits[1,2],  limits[2,1] <= y <= limits[2,2]
#  Plot coordiates are in pixels with 0 <= X <= X_max,  0 <= Y <= Y_max
function plot2data(X,Y,lims,X_max,Y_max)
	# convert plot coordinates (X,Y) to data coordinates (x,y)
	x = ((lims[1,2] - lims[1,1])/X_max) .*X .+ lims[1,1]
	y = ((lims[2,2] - lims[2,1])/Y_max) .*Y .+ lims[2,1]
	return (x,y)
end

function grid!(p_1,direc, limits, X_max, Y_max, col)
	# This function draws the background grid lines.
	# We first convert data coordinates to plotting coordinates to get grid lines
	# that are nicely spaced to the eye, regardless of actual data scaling.
	Δ = max(X_max,Y_max)/8
	Xlimits = [0,X_max]
    for (i,θ) in enumerate(direc)
		# Here θ is from data coordinates and Φ is from plotting coordinates. 
		# We work in plot coords to determine the grid lines, then convert to data
		# coords to do the actual plotting.
		tanϕ = tan(θ) * ((limits[1,2] - limits[1,1]) / (limits[2,2] - limits[2,1]) * Y_max / X_max)
		ϕ = atan(tanϕ)
		dX = Δ / abs(sin(ϕ))
		dY = Δ / abs(cos(ϕ))
        if tanϕ < 0
			lines = gridlines([0.0,X_max],Y_max,tanϕ,dX,dY)
			for line in lines
				(x,y) = plot2data(line[1],line[2],limits,X_max,Y_max)
				plot!(p_1,x,y,color=col[i],alpha=0.2)
			end
		else
			# we reflect the rectangle in X-direction and call gridlines, the reflect back
			lines = gridlines([-X_max,0.0],Y_max,-tanϕ,dX,dY)
			for line in lines
				(x,y) = plot2data(-line[1],line[2],limits,X_max,Y_max)
				plot!(p_1,x,y,color=col[i],alpha=0.2)
			end
		end
    end
end

function gridlines(Xlims,Ymax,tanϕ,dX,dY)
	# This function takes a rectangle defined by 
	#              Xlims[1] <= X <= Xlims[2],   0 <= Y <= Ymax
	# and draws grid lines dX or dY apart from each other along the coordinate axes. It
	# starts in the lower left corner and moves along the bottom then up the right side.  It
	# assumes ϕ is in [π/2,π), so tanϕ is negative and the grid lines are moving up and to 
	# the left.  The grid lines are defined as a vector of tuples (x,y), where x and y are 
	# 2-element vectors representing points on the boundary of the rectangle which connected 
	# make the grid line.
	lines = [] 
	xstart = min(Xlims[1] + dX,Xlims[2])
	ystart = 0.0
	# move along bottom to the right, then along right side upward
	while ystart < Ymax
		temp = xstart + (Ymax - ystart)/tanϕ
		if temp < Xlims[1]
			xend = Xlims[1]
			yend = ystart + tanϕ*(xend - xstart)
		else
			xend = temp
			yend = Ymax
		end
		push!(lines,([xstart,xend],[ystart,yend]))
		if xstart < Xlims[2]
			xstart += dX
			if xstart >= Xlims[2]  # end of bottom, now go up right side
				ystart = tanϕ * (Xlims[2] - xstart)
				xstart = Xlims[2]
			end
		else
			ystart += dY
		end
	end
	return lines
end

function iso_pre!(p_1,FG, direc, x, y, col)
	# plot the isoclines and their preimages for each θ in direc.
    (val_F,val_G) = FG(x,y)
    (val_FF,val_GG) = FG(val_F,val_G)
    for (i,θ) in enumerate(direc)
        # plot isoclines
        C = (val_F .- x) .* sin(θ) .- (val_G .- y) .* cos(θ) 
        iso = Contour.contour(x, y', C, 0.0)
        for line in lines(iso)
            xs, ys = coordinates(line) 
            plot!(p_1, xs, ys, linestyle=:solid, lw=1,color=col[i])
        end 

        #plot preimages of isoclines
        C = (val_FF .- val_F) .* sin(θ) - (val_GG .- val_G) .* cos(θ)
        pre = Contour.contour(x, y', C, 0.0)
        for line in lines(pre)
            xs, ys = coordinates(line)
            plot!(p_1, xs, ys, linestyle=:dash, lw=1,color=col[i])
        end

    end
end

function iso_pre_edge!(p_1, z, C_edge, LR_edge, val, linestyle, col)
	#This functon is to draw the isoclines and preimages on the edges.
	#iso_pre function fails to draw these because the contour package we use sometimes does not always
    #identify isoclines right along the edge.
	start=0
	stop=0
	#iterating through edge to find 0 values. start==0 means it has not started yet. 
	for i in 1:length(C_edge)
		if start == 0 && abs(C_edge[i]) < 1.0e-10
			start = i 
		elseif start > 0 && abs(C_edge[i]) >= 1.0e-10
			stop = i-1
		elseif start > 0 && length(C_edge) == i
			stop = i
		end
		if stop >start
			if LR_edge
				plot!(p_1,[val,val],[z[start],z[stop]],linestyle=linestyle, lw=1,color=col)
			else 
				plot!(p_1,[z[start],z[stop]],[val,val],linestyle=linestyle, lw=1,color=col)
			end
			start = 0
			stop = 0
		end
	end
end

function forwardmap!(p_1,FG,x,y) 
    (val_F,val_G) = FG(x,y)
    scatter!(p_1,val_F[:], val_G[:], color = :gray, marker = :circle, 
			    markersize = 1, alpha=0.008)
end

function detplotter!(p_1,x,y,FG,FG_Jac,col,x2,y2)
	# Points in (x,y) where det(J) is zero are plotted.
	# The image of points in (x2,y2) where det(J)=0 are plotted.
    det_J = Matrix{Float64}(undef,length(x),length(y))
    for j in 1:length(y), i in 1:length(x)
		det_J[i,j] = det(FG_Jac(x[i],y[j]))
    end
    det_zero = Contour.contour(x, y', det_J, 0)
    for line in lines(det_zero)
		 #Plot the determinant of the Jacobian =0 contour
		 xs, ys = coordinates(line)
		 plot!(p_1,xs, ys, linestyle=:dash, lw=1, color=col)
    end    
	if !isnothing(x2)
		det_J = Matrix{Float64}(undef,length(x2),length(y2))
		for j in 1:length(y2), i in 1:length(x2)
			det_J[i,j] = det(FG_Jac(x2[i],y2[j]))
		end
		det_zero = Contour.contour(x2, y2', det_J, 0)
	end
    for line in lines(det_zero)
		xs, ys = coordinates(line)
		#Plot the image of the determinat of the Jacobian =0 contour
		(val_F, val_G) = FG(xs, ys)
		plot!(p_1,val_F, val_G, linestyle=:solid, lw=1, color=col)
    end
end

function plot_annotate!(p_1,limits,FG,direc,labelpoints)
    # Analyze signs of regions in the main plot
	if !isnothing(labelpoints)
		points=labelpoints
	else
	    deltax = (limits[1,2] - limits[1,1])/10
	    deltay = (limits[2,2] - limits[2,1])/10
	    points = [(limits[1,1] + deltax, limits[2,1] + deltay), 
			      (limits[1,1] + deltax, limits[2,2] - deltay), 
			      (limits[1,2] - deltax, limits[2,1] + deltay), 
			      (limits[1,2] - deltax, limits[2,2] - deltay)]	
	end
    left_labels = label(FG, direc, points, iso=1)
    right_labels = label(FG, direc, points, iso=2)
	font_size = 9
    for (i,pt) in enumerate(points)
        (xi,yi) = pt
        #ploting "label_1 → label_2"
		scatter!(p_1,[xi],[yi],markersize=2,color=:black)
        annotate!(p_1,xi, yi, text("$(left_labels[i]) → $(right_labels[i])", :bottom, 
								   :hcenter, font_size, :black))
    end
end

function label(FG, direc, points; iso)
	# Determine the appropriate label for this point
    labels = Vector{String}(undef, length(points))
    
    for (j, (xi, yi)) in enumerate(points)
        dot_product = Vector{Float64}(undef,2)
        for (i,θ) in enumerate(direc)
            if iso==0
                dot_product[i] = xi .* (-sin(θ)) + yi .* cos(θ)
            elseif iso == 1
				(val_F,val_G) = FG(xi,yi)
                dot_product[i] = (val_F .- xi) .* (-sin(θ)) + (val_G .- yi) .* cos(θ)
            else  # iso=2
				(val_F,val_G) = FG(xi,yi)
				(val_FF,val_GG) = FG(val_F,val_G)
                dot_product[i] = (val_FF .- val_F) .* (-sin(θ)) + (val_GG .- val_G) .* cos(θ)
            end
        end
        if dot_product[1] >= 0 && dot_product[2] <= 0
            labels[j] = "1"
        elseif dot_product[1] >= 0 && dot_product[2] > 0
            labels[j] = "2"
        elseif dot_product[1] < 0 && dot_product[2] > 0
            labels[j] = "3"
        else  # dot_product[1] < 0 && dot_product[2] <= 0
            labels[j] = "4"
        end
        
    end
    return labels
end

function legendplot(limits,X_max,Y_max,direc,FG,col)
	# plot the legend
	# The legend data limits will be ±δx and ±δy.
	# The plot (pixel) limits will be ±ΔX and ±ΔY
	δx = (limits[1,2] - limits[1,1])/2.0
	δy = (limits[2,2] - limits[2,1])/2.0
	ΔX = X_max/2.0
	ΔY = Y_max/2.0
	# S is the factor for converting angles in data coordinates to plot coordinates
	S = ΔY*δx/(ΔX*δy)
	# A are the angles direc, converted to plot coordinates
	A = atan.(S.*tan.(direc))
	u = (A[1] + A[2])/2.0
	# phi are the angles in plot coords at which we want to plot the labels.
	ϕ = [u, u+π/2, u+π, u+3*π/2]
	
    # set up the plot
    p_2 = plot(xlims=(-δx, δx), ylims=(-δy, δy), framestyle=:none, axes=false, 
			   grid=false, legend=false)
    # lines of the cross which make θ angle with +x axis
    legendgrid!(p_2,direc, δx, col)

    # Add numbers in the space between lines
    font_size = 8
	r = 0.6
	plot_positions = [r*ΔX * cos.(ϕ)  r*ΔX * sin.(ϕ)]
	data_positions = eachrow([plot_positions[:,1].*(δx/ΔX)  plot_positions[:,2].*(δy/ΔY)])
	numbers = label(FG, direc, data_positions, iso=0)
    for (i,position) in enumerate(data_positions)
        annotate!(p_2,position[1], position[2], text(numbers[i], :black, :center, font_size))
    end
	return p_2
end

function legendgrid!(p_2,direc, x_1, col)
	# Draw the grid lines in the legend plot
    for (i,θ) in enumerate(direc)
        y_1 = x_1 * tan(θ)
        plot!(p_2,[-x_1, x_1], [-y_1, y_1], color=col[i], alpha=0.2)
    end
end
        
end
