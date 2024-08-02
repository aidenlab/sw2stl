# sw2stl

### The  project focuses on developing a Python script to convert binary files of virtual genomic structures from the “Spacewalk” web browser into STL files for 3D printing (or any other desired usage).

The Spacewalk browser, developed by the Aiden Lab, is a tool that uses igv.js and juicebox.js instances for rapid visualization and interaction between 3D and 1D genomic data.
Spacewalk's binary file format (.swb), later converted into .sw, is based on Hierarchical Data Format (HDF5), storing spatial data for genomic regions in the form of arrays: 
  - The genomic position group stores data such as chromosome names and start/end coordinates.
  - The spatial position group, crucial for this project, contains XYZ coordinate arrays linked to the genomic regions.

The program, in simple terms, follows these steps:
  - Parses .json file and extracts the URL of the .sw file from Spacewalk. It is downloaded if not available locally
  - Reads and extracts spatial position group XYZ coordinates from the .sw file (HDF5)
  - Spatial data is filtered to remove NaN values and XYZ coordinates are normalized to correct scale and position
  - Spheres are created at each coordinate point and cylinders connecting each consecutive point
  - All spheres and cylinders are combined into a single mesh. The mesh is checked for holes or self-intersections
  - After combination and validation, critical for successful 3D printing, the mesh is exported as an .stl file
  
  Why STL? Because it is easy to use, handles high complexity structures, and is a universal file format for everything 3D modeling.

  Future plans include refinement of algorithm efficiency, improvement of model quality and printing ease, exploration of multi-color models, embeding the program in Spacewalk, and development of additional conversion format options for users.
  