#%%
import argparse
import h5py
import numpy as np
import json
import wget
import re
import urllib.request
import trimesh
import os
from tqdm import tqdm

def get_sw_json(json_path):
    """Extract .sw URL from the .json file"""
    with open(json_path) as json_file:
        sw_json = json.load(json_file)
        sw_url = sw_json['spacewalk']['url']
    return sw_url

def download_sw_file(local_filename, sw_url):
    """Download the .sw file from the SpaceWalk URL"""
    urllib.request.urlretrieve(sw_url, local_filename)

def create_sphere(center, radius):
    """Create a sphere at the given center with the given radius"""
    sphere = trimesh.creation.icosphere(subdivisions=2, radius=radius) 
    sphere.apply_translation(center)
    return sphere

def create_cylinder(start_point, end_point, radius):
    """Create a cylinder with the given start and end points and radius"""
    height = np.linalg.norm(np.array(end_point) - np.array(start_point))
    direction = np.array(end_point) - np.array(start_point)
    direction /= height
    
    cylinder = trimesh.creation.cylinder(radius=radius, height=height, sections=16) 
    
    # Calculate rotation matrix
    z_axis = np.array([0, 0, 1])
    v = np.cross(z_axis, direction)
    c = np.dot(z_axis, direction)
    k = 1 / (1 + c)
    
    rotation = np.array([
        [v[0] * v[0] * k + c, v[0] * v[1] * k - v[2], v[0] * v[2] * k + v[1]],
        [v[1] * v[0] * k + v[2], v[1] * v[1] * k + c, v[1] * v[2] * k - v[0]],
        [v[2] * v[0] * k - v[1], v[2] * v[1] * k + v[0], v[2] * v[2] * k + c]
    ])
    
    # Convert rotation matrix to 4x4 transformation matrix
    transformation = np.eye(4)
    transformation[:3, :3] = rotation
    cylinder.apply_transform(transformation)
    cylinder.apply_translation((start_point + end_point) / 2.0)
    
    return cylinder

def normalize_vertices(vertices):
    """Normalize vertices to fit within a unit cube centered at the origin"""
    min_coords = np.min(vertices, axis=0)
    max_coords = np.max(vertices, axis=0)
    center = (min_coords + max_coords) / 2
    vertices -= center  # Center the vertices
    scale = np.max(max_coords - min_coords)  # Find the scale based on the largest dimension
    vertices /= scale  # Normalize to fit within a unit cube
    return vertices, scale

def convert_to_stl_and_offset(hdf5_file_path, ensemble_group, output_stl, time_point):
    """Convert the spatial position data from the hdf5 file to an .stl file with spheres and cylinders"""
    # Open the HDF5 file and extract the spatial position data (vertex coordinates)
    with h5py.File(hdf5_file_path, 'r') as swbf:
        dataset = swbf[ensemble_group]
        data = dataset['spatial_position']
        sp = data[time_point]
        vertices = np.array(sp, dtype=np.float32)

    # Process and refine vertices (remove NaNs)
    vertices = vertices[~np.isnan(vertices).any(axis=1)]

    # Normalize the vertices
    vertices, scale = normalize_vertices(vertices)

    # Set dynamic sphere and cylinder sizes based on the normalized scale
    sphere_radius = 0.01  # Small radius for spheres
    cylinder_radius = 0.005  # Small radius for cylinders

    # Initialize an empty mesh
    combined_mesh = None

    for i in tqdm(range(len(vertices))):
        sphere = create_sphere(vertices[i], sphere_radius)
        
        if combined_mesh is None:
            combined_mesh = sphere
        else:
            combined_mesh = trimesh.boolean.union([combined_mesh, sphere])

        if i < len(vertices) - 1:
            cylinder = create_cylinder(vertices[i], vertices[i + 1], cylinder_radius)
            combined_mesh = trimesh.boolean.union([combined_mesh, cylinder])

    # Ensure the combined mesh is watertight (no holes or self-intersections)
    if not combined_mesh.is_watertight:
        combined_mesh = combined_mesh.fill_holes()

    # Save the final combined mesh as an .stl file
    combined_mesh.export(output_stl, file_type='stl')
    print(f"STL file saved as {output_stl}")

def main():
    """Define arguments, determine .json or .sw to proceed, and call the convert_to_stl_and_offset function"""
    parser = argparse.ArgumentParser(description='Convert dataset to STL and apply offset.')
    parser.add_argument('-f','--file', required=True, help='Path to the file')
    parser.add_argument('-j','--json_file', required=False, help='Path to the json file')
    parser.add_argument('-hd','--hdf5_file', required=False, help='Path to the hdf5 file')
    parser.add_argument('-e','--ensemble_group', required=True, help='Title of the ensemble group (cap senstive)')
    parser.add_argument('-t','--time_point', required=True, help='''Time point number. 
                        Note: The time point should be in the format "t_#" where # is the time point number. 
                        Also consider model 1 = t_0, therefore, may need to -1 the # before plugging in to ensure correct model''')
    parser.add_argument('-o','--output_stl', default=None, help='Output STL file title (Optional)')
    
    args = parser.parse_args()
    
    hdf5_file_path = None
    
    if '.json' in args.file:
        args.json_file = args.file
        
        sw_url = get_sw_json(args.json_file)
    
        if 'dropbox' in sw_url:
            if 'dl=0' in sw_url:
                sw_url = sw_url.replace('dl=0', 'dl=1')
            hdf5_file_path = wget.download(sw_url)
        else:
            pattern = r'/([^/]+)$'
            match = re.search(pattern, sw_url)
            if match:
                extracted_part_url = match.group(1)
            else:
                raise ValueError("No match found in the URL")
            extracted_part_url = extracted_part_url.replace('.sw', '')
            hdf5_file_path = f'{extracted_part_url}.sw'
            if hdf5_file_path in os.listdir():
                print(f"{hdf5_file_path} already exists in the current directory")
            else:
                download_sw_file(hdf5_file_path, sw_url)
            
    elif '.sw' in args.file:
        args.hdf5_file = args.file
        hdf5_file_path = args.hdf5_file
    else:
        raise ValueError("Please provide the path a file that is either a .json or .sw")
    
    if args.output_stl is None:
        args.output_stl = f'{args.ensemble_group}_{args.time_point}.stl'

    
    convert_to_stl_and_offset(hdf5_file_path, args.ensemble_group, args.output_stl, args.time_point)

if __name__ == '__main__':
    main()



# example parameters

    # Mammoth: -j = ./mammoth-spacewalk-session.json, -e = replica01_chr10_Direct_Inversion, -t = t_141
    # HCT116 6h_auxin chr21:28-30: -j = ./HCT116-spacewalk-session.json, -e = HCT116, -t = t_14
    # HCT116 6h_auxin chr21:28-30: -j = ./HCT116-t-36-spacewalk-session.json, -e = HCT116, -t = t_36
    # K562 chr21:28-30: -j = ./K562-t-1008-spacewalk-session.json, -e = K562, -t = t_1008


