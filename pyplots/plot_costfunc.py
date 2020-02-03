from mpl_toolkits import mplot3d # can now pass kw projection="3d" to any of the normal axes creation routines
import mpl_toolkits.mplot3d.axes3d as p3 # for contourf3D
import matplotlib.pyplot as plt
import numpy as np

# an aesthetic function
aestheticf = lambda x: np.sin(x[0])**10 + np.cos(10.+x[1]*x[0])*np.cos(x[0])
aestheticfdom = ((0.,5.),(0.,5.))

# six-hump camel function (x is 2D args)
shcamelf = lambda x: (4. - (2.1*x[0]**2) + ((x[0]**4)/3.))*x[0]**2 + (x[0]*x[1]) + \
                    (-4. + (4.*x[1]**2))*x[1]**2

# PARAMETERS
ndim = 2 # dimensionality of cost function
costfunc = shcamelf # handle for cost function to be plotted
domain = ((-2.,2.),(-1.,1.)) # domain range ((x1_min,x1_max),(x2_min,x2_max),...)
'''
costfunc = aestheticf
domain = aestheticfdom
'''
#resln = 0.01
resln = 0.1

assert len(domain) == ndim

d2=False
d3=True

# EXAMPLE PLOTTING OF 2D FUNCTION AS CONTOUR PLOT
if d2:
    print "Plotting 2d contour plot of function"
    plt.style.use("seaborn-white") # style

    x_no_samples = int((domain[0][1]-domain[0][0])/resln)
    y_no_samples = int((domain[1][1]-domain[1][0])/resln)
    x = np.linspace(domain[0][0],domain[0][1],x_no_samples)
    y = np.linspace(domain[1][0],domain[1][1],y_no_samples)
    X, Y = np.meshgrid(x,y)
    Z = costfunc((X,Y))
    '''
    # option one - unfilled colour
    contours = plt.contour(X,Y,Z,6,colors="black") # unfilled colour
    plt.clabel(contours,inline=True,fontsize=8) # label the contours
    plt.imshow(Z,extent=[domain[0][0],domain[0][1],domain[1][0],domain[1][1]],origin="lower",
               cmap="RdGy",alpha=0.5) # alpha is transparency
    '''
    # option two - filled colour
    plt.contourf(X,Y,Z,50,cmap="RdGy") # filled colour

    plt.colorbar()
    plt.show()


# EXAMPLE PLOTTING OF 2D FUNCTION IN THREE DIMENSIONS
elif d3:
    print "Plotting 3d projection of function"
    x_no_samples = int((domain[0][1]-domain[0][0])/resln)
    y_no_samples = int((domain[1][1]-domain[1][0])/resln)


    x = np.linspace(domain[0][0],domain[0][1],x_no_samples)
    y = np.linspace(domain[1][0],domain[1][1],y_no_samples)
    X, Y = np.meshgrid(x,y)
    Z = costfunc((X,Y))

    '''
    # option one - 3d contour plot
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.contour3D(X,Y,Z,50,cmap="coolwarm")
    ax.set_xlabel("x1")
    ax.set_ylabel("x2")
    ax.set_zlabel("V(x1,x2)")
    ax.set_title("cost function")
    ax.view_init(60,35) # set elevation and azimuthal angles for optimal viewing
    '''
    '''
    # option two - wireframe plot (set resolution to low e.g. 0.2)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_wireframe(X, Y, Z, color='black')
    ax.set_title('wireframe');
    '''
    #'''
    # option three - surface (like filled wireframe) plot
    # set resolution to low (e.g. 0.1 or 0.2)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    # note ordering: higher zorder means on top
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none',zorder=1)
    # hide grid and axes
    plt.axis("off")
    ax.grid(False)
    #'''
    #'''
    # OTHER STUFF - ADD LINES AND POINTS
    # plot random points
    n_points = 15
    np.random.seed(19)
    xdata = np.random.uniform(domain[0][0],domain[0][1],n_points)
    ydata = np.random.uniform(domain[1][0],domain[1][1],n_points)
    zdata = np.array([costfunc((x1,x2)) for (x1,x2) in np.column_stack((xdata,ydata))])
    ax.scatter3D(xdata,ydata,zdata,c=zdata,cmap="Reds",marker="^",s=200,zorder=2)
    # plot a line on the surface (linear interpolation) connecting the two global minima
    min1, min2 = (0.0898,-0.7126), (-0.0898, 0.7126)
    xline = np.linspace(min1[0],min2[0],1000)
    yline = np.linspace(min1[1],min2[1],1000)
    zline = np.array([costfunc((x1,x2)) for (x1,x2) in np.column_stack((xline,yline))])
    ax.plot3D(xline,yline,zline,"gray",linewidth=2.,zorder=3)
    #'''
    '''
    # option four - surface 3D plot with partial polar (not rectilinear) grid
    r = np.linspace(0, 6, 20)
    theta = np.linspace(-0.9 * np.pi, 0.8 * np.pi, 40) # partial polar - give 'slice' of function
    # theta = np.linspace(-np.pi,np.pi,40) # alternative: polar - view complete function
    r, theta = np.meshgrid(r, theta)

    X = r * np.sin(theta)
    Y = r * np.cos(theta)
    Z = costfunc((X, Y))

    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none');

    '''
    '''
    # option five - filled 3D contour plot
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.contourf3D(X,Y,Z,50)
    fig.add_axes(ax)
    '''

    plt.show()
