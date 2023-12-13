# Instruction to create a CAD design and a mesh

## Tetraedra

### CAD design
One of the open source option is [FreeCAD](https://www.freecadweb.org/). To install it, follow the instructions given [here](https://www.freecadweb.org/wiki/Installing).  

To learn how to use it, [this](https://www.youtube.com/channel/UCIJuGbJfjFqCIPcEi2EbVHw) youtube channel is usefull. 
Moreover, there are plenty of examples and other tutorials that can be found online.  

Once you have a CAD design, save it using the .brep extension (go under file -> export -> and choose .brep). A toy example can be found in the directory in the folder example/.

### Mesh
To do the tetraedra meshing, [gmsh](http://gmsh.info/) is a good option.  

It can be downloaded following the instructions [here](http://gmsh.info/#Download).  

To do a mesh, the easiest way is to use the terminal and type the following command:

```
gmsh path/to/file -int_dimension -o name -format extension -clmin min_len -clmax max_len
```
Where *int_dimension* is the number of dimensions, *extension* the extension of the output file, *min_len* the minimum mesh element size and *max_len* the maximum mesh element size.  

For example, if you run the command in the directory where you have the example CAD file sphere.brep, and you want to do a 3D mesh 
with both the minimum and maximum mesh size to be 20, an output format .msh and the name sphere, you would run the following command:

```
gmsh sphere.brep -3 -o sphere -format msh -clmin 20 -clmax 20
```
To visualize the meshing, you can use gmsh by running the following command:

```
gmsh name
```
In our example:

```
gmsh sphere
```

If needed, all the command line options can be found [here](http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options).  

Also, it is important to pay attention to the written output in the terminal after the meshing process. It contains usefull information about the meshing. Notably, if for some geometrical reasons the meshing cannot be done, it will be written there.  

Finally, you can also open the .msh file with your favorite text editor to see the details. The instructions to read this file can be found [here](http://gmsh.info//doc/texinfo/gmsh.html#MSH-ASCII-file-format).
