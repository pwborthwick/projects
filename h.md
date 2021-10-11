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

    parameters - none. This generates random values of l, m, &theta; and &phi;, and then calculates the (scipy) value from angularSolution routine comparing it with the (sympy) value from the **Ynm** function. Returns True or False.

7.  **radialSolution(n, l, r)**

    parameters - *n* is the principal quantum number, *l* is the angular quantum number (0 &ge; l &le; n) and *r* is the radial distance (in units of the Bohr radius). Calls scipy.special.genlaguerre to return the normalised radial wavefunction value at r.

8.  **radialVerify()**

    parameters - none. This generates random values of n, l and r , and then calculates the (scipy) value from radialSolution routine comparing it with the (sympy) value from the **R_nl** function. Returns True or False.

9.  **radialSolutionType(n, l, r, psi_type = 'radial distribution')**

    parameters - *n* is the principal quantum number, *l* is the angular quantum number (0 &ge; l &le; n), *r* is the radial distance (in units of the Bohr radius) and *psi_type* is the type of wavefunction to be plotted. *psi_type* can be one of 'radial distribution' (&psi;), 'probability density' (|&psi;|<sup>2</sup>) and 'probability distribution' (4&pi;r<sup>2</sup>|&psi;|<sup>2</sup>). Returns the value of the *psi_type*.

10. **radialSolutionPlot(n, l, psi_type = 'radial distribution', psi_normal = False, parameters = {'points':100, 'size':[7,5], 'extent':[20,0.2], 'equal':False})**

    parameters - *n* is a list of principal quantum numbers to plot, *l* is a list of angular quantum numbers corresponding to the list of principal quantum numbers (0 &ge; l &le; n) to plot, *r* is the radial distance (in units of the Bohr radius) and *psi_type* is the type of wavefunction to be plotted. *psi_type* can be one of 'radial distribution' (&psi;), 'probability density' (|&psi;|<sup>2</sup>) and 'probability distribution' (4&pi;r<sup>2</sup>|&psi;|<sup>2</sup>). *psi_normal* determines if the wavefunction should be normalised. *parameters* is a dictionary of values defining the plot appearance, these values are 'points' - the number of data points to calculate, 'size' - the width and height of the plot in inches, 'extent' - the horizontal and vertical ranges and 'equal' - forces a square aspect ratio. Plots the radial distribution. This is the plot for\
    radialSolutionPlot([1,2,3,2,3,3], [0,0,0,1,1,2], 'radial distribution', False, {'points':100,'size':[7,5],'extent':[20,[-0.15,0.2]],'equal':False})
    
![image](https://user-images.githubusercontent.com/73105740/136762722-b95f78b7-3f3a-4a40-9898-f2aa2119c2c1.png)

11. **wavefunction(n, l, m, r, theta, phi, grid)**

    parameters - *n* is a principal quantum numbers to plot, *l* is an angular quantum numbers, *m* is the magnetic quantum number, *r*, *theta* and *phi* are the coordinates of a point in spherical polars or cartesian (*grid*='cartesian'). Calculates the product of radial and angular solutions and returns value.

12. **wavefunctionContour(n, l, m, parameters={'points':100, 'extent':[-20, 20], 'color_map':'coolwarm','plane':'xy', 'elevation':0.0, 'contour': False})**

    parameters - *n* is a principal quantum numbers to plot, *l* is an angular quantum numbers, *m* is the magnetic quantum number and *parameters* is a dictionary of values defining the plot appearance, these values are 'points' - the number of data points to calculate, 'extent' - the horizontal and vertical ranges, 'color_map' - is the matplotlib name of the color theme to use, 'plane' - specifies which plane to view the contour in, 'elevation - is the height of the slice in 'plane' and 'contour' - specifies if the contour lines should be plotted. Shows a contour plot of the wavefunction. Examples for\
    wavefunctionContour(3, 2, 0, {'points':80, 'extent':[-30, 30], 'color_map':'gist_yarg', 'plane':'zx', 'elevation':0, 'contour': True})\
	wavefunctionContour(3, 2, 0, {'points':80, 'extent':[-30, 30], 'color_map':'gray', 'plane':'zx', 'elevation':0, 'contour': False})
![image](https://user-images.githubusercontent.com/73105740/136765639-5d801e98-995e-4d87-8760-6ba1925fe9a2.png)

13. **wavefunctionVerify()**

    parameters - None. This generates random values of *n*, *l*, *m* and *r*, *theta*, *phi* and uses these to compare wavefunction values generated by the routine wavefuction and the (sympy) routine Psi_nlm. Returns True or False.


