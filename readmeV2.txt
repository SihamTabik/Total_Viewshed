This code has been created for total viewshed computation, but the same kernel, with small changes,
can also be used for related computations, including Total-3D-Viewshed, single point (tower) or many
point (track) viewshed, etc.

The model can only deal (now) with square DEMS, which must be provided by either 16-bit .bil (short) 
or 32-bit .flt (float) format files. A .hdr file with nrows=ncols and cellsize is recommended
although the first parameters should be detected from the file size (no gaps in the files considered).

By default, observers are placed 1.5 (meters) above the DEM (the position of the programmer's eyes :).

A boolean mask can also be provided to determine an area-of-interest. Note that points with mask=true
are considered of null interest. 


Usage: tvs.exe model [-v] [-q] [-ttype][-Hheight][-mmaskfilename] [-nitem]
            -v Verbose mode
            -q Silent mode
            -tN Execution mode:
              N=0 Total viewshed (default)
              N=1 Total 3D-Viewshed
              N=2 Single tower viewshed
              N=3 Track Viewshed
              N=4 Single tower 3D-viewshed
              N=5 Track 3D-Viewshed
              N=6 Find sequential towers (cell coverage algorithm)
              N=7 Isolated area detection
              N=8 Horizon computation
              N=9 Timing mode
            -Hheight Observer's height in meters (default, 1.5m)\n"
            -m Use default mask file (input/mask.dat for area of interest\n"
            -mfile Use file for area of interest\n"
            -nitem In single-tower and several-tower (tracks) viewshed execution modes, specifies the track:\n"
              item:utmn,utme,name\n"
              item:index,name\n"
              item:index\n"
              Use custom function to assign indexes to towers\n\n"

Example:

 ./comp.out 4075000_0310000_010.bil -v -minput/mask.dat -t2 -n4063520,312340,mytower1

model:   The filename .bil or .flt. without route. These model must be placed in a directory ./input. Some sample DEMs 
         are provided in little endian format, so the code should only work with these DEMs in little endian 
         (intel) machines. A .hdr file with nrows and ncols is recommended.
         
         Model data with a single band (elevation) in .bil format are considered only as unsigned short integers in 
         natural ordering (the inner loop runs from east to west, while the outer loop goes from north to south. In .flt 
         files, the same ordering is used, but elevation data are IEEE float numbers.



-v:      Exhaustive information.

-q:      Quiet mode, without messages

-tN:     The type of algorithm. Two groups of algorithms are considered:

         Total algorithms (type 0 and type 1). These algorithms calculates 2d and 3d viewshed for all point dems,
         and takes full profit of our algorithm. With a 12 core xeon architecture, it requires about five minutes
         for 4000000 potential observers. There is a flag (fulldata) which is used to store the geometry of all
         viewsheds, but this is usually unnecessary.
         
         Single point or multiple point viewshed. In single point 2D-viewshed, the algorithm is probably similar or
         worse than other methods (GRASS, ArcGIS), and take a few seconds to compute the results. Otherwise, for 
         thousands of points, we use the same kernel that total viewsheds and you can take the results in a few minutes.
         
         In the second group of algorithm, you can provide you own set by changing some routines in auxfunc.cpp 
         The routines uses the index provide with the next flag (-n) for the selection of candidate observers.
         The sample routines uses either a set of points given in a text file (track.txt) or a binary file (list.dat)
         
-nitem   In single-tower and several-tower (tracks) viewshed execution modes, specifies the point or set.
         For a single point, the format can be: -nutmn,utme,name
         Two other format are -nindex and -nindex,name. In this cases, you should change to customize selection 
         functions (setCustomTowers, setCustomTracks)

-Hh:     Observer's height h in default unit (meters) (default, 1.5m)

-mfile: Mask file. The size (nros, ncols) must match the size of the DEM. Natural ordering is used as described above.

-?       Help

Some sample models are provided:

4070000_0310000_010.bil (2000x2000) Parque Natural Sierra de las Nieves
4080000_0360000_010.bil (2000x2000) Municipality of Malaga
4075000_0310000_010.bil (2500x2500) Extended area Sierra de las Nieves


