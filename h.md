## Hydrogenic Wavefunctions

1.  **transform(mode, i, j, k)**

    parameters - *mode* is either 'cartesian->spherical' or 'spherical->cartesian', *i* is either an x-coordinate or an r-coordinate, *j* is either a y-coordinate or a &theta;-coordinate and *k* is either a z-coordinate or a &phi;-coordinate depending on the *mode*. The routine performs a coordinate transformation according to the *mode* and returns the transformed coordinates.

2.  **angularSolution(m, l, theta, phi)**

    parameters - *m* is the magnetic quantum number, *l* is the angular quantum number, *theta* is the polar coordinate and *phi* the azimuthal coordinate. Calls routine scipy.special.sph_harm and returns the spherical harmonic Y<sup>m</sup><sub>l</sub>(&theta;,&phi;).

3.  **angularSolutionPlot(ax, m, l, parameters)**

    parameters - *ax* are the matplotlib axes as define eg from 'fig.add_subplot', *m* is the magnetic quantum number, *l* is the angular quantum number, and *parameters* is a directory of values defining the plot. The *parameters* directory is defined as follows - parameters = {'points':70,'extent':[-0.5,0.5],'color_map':'coolwarm', 'bar':'on','axes':'off','alpha':0.8} where 'points' are the number of data points in the (&theta;, &phi;)- grid, 'extent' is extent of the radial values to display, 'color_map' is the matplotlib color mapping to use, 'bar' determines whether to draw a heat bar, 'axes' controls whether to draw the coordinate axes or not and 'alpha' is the transparancy factor. Uses matplotlib.plot_surface to display the spherical harmonic in the color map requested. Optionally displays a colorbar. Does not 'show' the plot.

4.  **angularSolutionPlotSingle(m, l, parameters = {'points':70,'extent':[-0.5,0.5],'color_map':'coolwarm','bar':'on','axes':'off','alpha':0.8})**

    parameters - *m* is the magnetic quantum number, *l* is the angular quantum number and and *parameters* is a directory of values defining the plot. This routine displays a spherical harmonic defined by *m* and *l* according to the display conditions in *parameters*. An example for Y</sub>3,0</sub> is
    
    ![image](https://user-images.githubusercontent.com/73105740/136229523-1f09d93e-3834-497a-a2e5-5bbbd8db491e.png)


5.  **angularSolutionPlotFamily(l_maximum, parameters = {'points':70,'extent':[-0.5,0.5],'color_map':'coolwarm','bar':'off','axes':'off','alpha':1.0})**

    parameters - *l_maximum* is the highest value of the angular quantum number to display. *parameters* is a directory of values defining the plot. This routine will display all spherical harmonics for l=0, l=1, ... l=l_maximum where -l &le; m &ge; l. An example for l_maximum = 3 is
    
    ![image](https://user-images.githubusercontent.com/73105740/136230126-100046a3-7ee5-4435-ae2f-dc43eb5c22b3.png)


6.  **angularVerify()**

    parameters - none. This generates random values of l, m, &theta; and &phi; and calculates the (scipy) value from angularSolution routine comparing it with the (sympy) value from Ynm function. Returns True or False.

    
