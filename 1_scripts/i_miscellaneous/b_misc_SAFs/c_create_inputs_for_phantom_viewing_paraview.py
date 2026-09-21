import meshio

# Read the TetGen mesh
mesh = meshio.read("MRCP_Female.ele")

# Save as a ParaView-compatible file
meshio.write("MRCP_Female.vtu", mesh)

print("Created MRCP_Female.vtu")