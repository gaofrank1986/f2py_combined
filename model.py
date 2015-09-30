import numpy as np
import te

te.mesh.read_mesh()
te.mesh.convsb()
te.mesh.prepare_mesh()
te.gen_matrix.calc_matrix()

print te.gen_matrix.amata.diagonal()
