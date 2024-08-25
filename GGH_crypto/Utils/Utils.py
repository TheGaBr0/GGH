from flint import fmpz_mat, fmpq_mat, fmpq, fmpz
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal, getcontext
import math
import os
import subprocess
import ast
import threading


class Utils:
    def gram_schmidt(matrix_fmpz):
        matrix = fmpq_mat(matrix_fmpz)
        n, m = matrix.nrows(), matrix.ncols()
        
        result = fmpq_mat(n, m)
        w_array = [[fmpq(0) for _ in range(m)] for _ in range(n)]

        for i in range(n):
            v_i = [matrix[i, j] for j in range(m)]
            w_i = v_i[:]
            
            for j in range(i):
                w_j = w_array[j]
                
                # Calculate dot products
                dot_v_w = sum(v_i[k] * w_j[k] for k in range(m))
                dot_w_w = sum(w_j[k] * w_j[k] for k in range(m))
                
                # Perform subtraction
                factor = dot_v_w / dot_w_w
                w_i = [w_i[k] - factor * w_j[k] for k in range(m)]
            
            w_array[i] = w_i
            
            for j in range(m):
                result[i, j] = w_i[j]

        return result
    
    def start_visualization(basis_1, basis_1_cvp, point, basis_2=None, basis_2_cvp=None, title="Lattice Plot", limit=5):
        thread = threading.Thread(target=Utils.generate_lattice_points, args=(basis_1, basis_1_cvp, point, basis_2, basis_2_cvp, title, limit))
        thread.start()
        return thread
    
    def generate_lattice_points(basis_1, basis_1_cvp, point, basis_2=None, basis_2_cvp=None, title="Lattice Plot", limit=5):

        basis_1_np = np.array(basis_1.tolist()).astype(int)

        # Check if basis_2 is provided
        if basis_2 is not None:
            basis_2_np = np.array(basis_2.tolist()).astype(int)

        point_np = np.array(point.tolist()).astype(int)

        # Create a meshgrid of integer coordinates
        x = np.arange(-limit, limit + 1)
        y = np.arange(-limit, limit + 1)
        coords = np.array(np.meshgrid(x, y)).T.reshape(-1, 2)

        # Multiply the meshgrid by the basis matrix to get lattice points
        lattice_points = np.dot(coords, basis_1_np)
        
        plt.plot(lattice_points[:, 0], lattice_points[:, 1], 'bo', color='black')
        plt.grid(True)
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.title(title)
        label_w = "w"
        label_w_2 = "w_2"
        label_t = "t"
        
        plt.scatter(point_np[0][0], point_np[0][1])
        
        if basis_2_cvp is not None:
            plt.scatter(int(basis_2_cvp[0, 0]), int(basis_2_cvp[0, 1]), color='red', s=70)
            
        if basis_1_cvp is not None:
            plt.scatter(int(basis_1_cvp[0, 0]), int(basis_1_cvp[0, 1]), color='blueviolet', s=70)
            
        plt.annotate(label_t, # this is the text
                    (point_np[0][0], point_np[0][1]), # these are the coordinates to position the label
                    textcoords="offset points", # how to position the text
                    xytext=(0, -12), # distance from text to points (x,y)
                    ha='center') # horizontal alignment can be left, right or center
                    
        if basis_1_cvp is not None:
            plt.annotate(label_w, 
                        (int(basis_1_cvp[0, 0]), int(basis_1_cvp[0, 1])), 
                        textcoords="offset points",
                        xytext=(0, 10), 
                        ha='center')
                        
        if basis_2_cvp is not None:
            plt.annotate(label_w_2, 
                        (int(basis_2_cvp[0, 0]), int(basis_2_cvp[0, 1])), 
                        textcoords="offset points", 
                        xytext=(0, -12), 
                        ha='center') 
            
        plt.arrow(0, 0, basis_1_np[0][0], basis_1_np[0][1], lw=2, head_width=0.1, head_length=0.1, fc='blueviolet', ec='blueviolet')
        plt.arrow(0, 0, basis_1_np[1][0], basis_1_np[1][1], lw=2, head_width=0.1, head_length=0.1, fc='blueviolet', ec='blueviolet')
        
        if basis_2 is not None:
            plt.arrow(0, 0, basis_2_np[0][0], basis_2_np[0][1], lw=2, head_width=0.1, head_length=0.1, fc='red', ec='red')
            plt.arrow(0, 0, basis_2_np[1][0], basis_2_np[1][1], lw=2, head_width=0.1, head_length=0.1, fc='red', ec='red')  
        
        plt.show()

    def embedding(basis, point, visualize=False, GGH=False):
        
        if basis.ncols() != point.ncols():
                raise ValueError(f"[Utils] Point is a {1}x{point.ncols()} matrix, but basis is a {basis.nrows()}x{basis.ncols()} one")
        
        n = basis.nrows()

        # Create the initial matrix with R and a column of zeros
        matrix_emb = fmpz_mat([[int(basis[i,j]) if j < n else 0 for j in range(n+1)] for i in range(n)])

        # Add t as the last row, with 1 as the last element
        last_row = [int(point[0,i]) for i in range(n)] + [1]
        matrix_emb = fmpz_mat(matrix_emb.tolist() + [last_row])
        
        matrix_emb = matrix_emb.lll()

        min_norm = float('inf')
        shortest_vector = None

        for i in range(n + 1):
            i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n + 1)]])

            if GGH:
                if matrix_emb[i, n] == 1:  # Se l'ultimo elemento è 1
                    norm = Utils.vector_norm(i_th_row)
                    if norm < min_norm:
                        min_norm = norm
                        # Prendi i primi n valori della riga più corta
                        shortest_vector = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])
            else:
                # Se GGH non è vero, considera comunque la riga per trovare il vettore più corto
                norm = Utils.vector_norm(i_th_row)
                if norm < min_norm:
                    min_norm = norm
                    # Prendi i primi n valori della riga più corta
                    shortest_vector = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])
    
        # Se GGH è vero e nessun vettore con 1 come ultimo elemento è stato trovato,
        # ritorna comunque il vettore più corto trovato
        if GGH and shortest_vector is None:
            for i in range(n + 1):
                i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n + 1)]])
                norm = Utils.vector_norm(i_th_row)
                if norm < min_norm:
                    min_norm = norm
                    # Prendi i primi n valori della riga più corta
                    shortest_vector = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])
  
        closest_vector = point - shortest_vector

        if visualize:
            if basis.nrows() != 2:
                raise ValueError(f"[Utils] Can't visualize. Basis must be a {2}x{2} matrix, but got a {basis.nrows()}x{basis.ncols()} one")
            thread = Utils.start_visualization(basis, closest_vector, point, title="Embedding method")
            thread.join()

        return closest_vector
        
    def vector_norm(row):
        if isinstance(row, fmpz_mat):
            row = fmpq_mat(row)
        getcontext().prec = 50
   
        return Decimal(sum((Decimal(int(x.numer())) / Decimal(int(x.denom()))) ** 2 for x in row)).sqrt()
    
    def get_hadamard_ratio(basis = None):
        matrix = basis
        norms = []
        dimension = basis.nrows()

        for i in range(matrix.nrows()):
            row = fmpz_mat([[matrix[i, j] for j in range(matrix.ncols())]])
            norm = Utils.vector_norm(row)
            norms.append(norm)
            
        denominator = math.prod(norms)
        numerator = abs(Decimal(matrix.det().str()))

        result = (numerator / denominator) ** Decimal(1 / dimension)
        return result, f"{result:.50f}"
    
    def write_matrix_to_file(matrix, filename):
        filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)

        rows = matrix.nrows()
        cols = matrix.ncols()

        # Open the file for writing
        with open(filename, "w") as file:
            file.write("[")
            # Iterate over each row of the matrix
            for i in range(rows):
                # Iterate over each column of the matrix
                file.write("[")
                for j in range(cols):
                    # Write the element to the file
                    file.write(matrix[i, j].str() + " ")
                # Write a newline character after each row
                file.write("]\n")
            file.write("]")
            
        return filename
    
    def load_matrix_from_file(filename, fplll = False):
        filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)

        if not fplll:

            with open(filename, 'r') as file:
                matrix_data = file.read().strip()

            # Removing the outer brackets
            matrix_data = matrix_data[1:-1]

            # Splitting into rows and then into individual elements
            rows = matrix_data.split('\n')
            matrix_list = [list(map(int, row.strip()[1:-1].split())) for row in rows]
            print(matrix_list)
            return fmpz_mat(matrix_list)
        else:
            with open(filename, 'r') as file:
                content = file.read()
            
            # Remove whitespace and newlines
            content = content.replace(']', '],')
            content = content.replace(' ]', ']')
            content = content.replace(' ', ', ')
            content = content[:-5] + ']'

            # Use ast.literal_eval to safely evaluate the string as a Python expression
            matrix = ast.literal_eval(content)
            
            return fmpz_mat(matrix)
    
    def BKZ(matrix, block=20, pruned=False):
        path = Utils.write_matrix_to_file(matrix, 'BE.txt')
        output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out.txt')

        if os.name == 'nt':
            path = '/mnt/' + path.replace(":", "")  # Windows path to WSL format
            path = path.replace("\\", "/")
            path = path.replace(" ", "\ ")

            output_path = '/mnt/' + output_path.replace(":", "")
            output_path = output_path.replace("\\", "/")
            output_path = output_path.replace(" ", "\ ")

            if not pruned:
                command = f"wsl fplll {path} -a bkz -b {block} -p 85 -f mpfr -m proved > out.txt"
            else:
                command = f"wsl fplll {path} -a bkz -b {block} -p 85 -f mpfr -s default.json -bkzmaxloops 3 > out.txt"
        else:
            command = f"fplll '{path}' -a bkz -b {block} -p 100 -f mpfr -m proved > out.txt"

        try:
            # Run the command and capture its output
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output, error = process.communicate()
            output = output.decode('utf-8')
            error = error.decode('utf-8')

            if error:
                print("Error during reduction:", error)

            # Move the output file to the desired location
            command = f"wsl mv -f out.txt {output_path}"
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output, error = process.communicate()
            error = error.decode('utf-8')

            if error:
                print("Error during file move:", error)

            print("Reduction completed")

            # Load the reduced matrix
            reduced_matrix = Utils.load_matrix_from_file("out.txt", True)

            # Remove the temporary files with error checking
            remove_command = f"rm -f {path} {output_path}" if os.name != 'nt' else f"wsl rm -f {path} {output_path}"
            process = subprocess.Popen(remove_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output, error = process.communicate()
            error = error.decode('utf-8')

            if error:
                print("Error during file removal:", error)

            return reduced_matrix

        except Exception as e:
            return None, str(e)
        
    def babai_rounding(basis, point, visualize=False):

        a = point * basis.inv()
       
        cols = a.ncols()
        rows = a.nrows()

        for i in range(rows):
            for j in range(cols):
                a[i,j] = round(a[i,j])
        
        closest_vector = a * basis 

        if visualize:
            if basis.nrows() != 2:
                raise ValueError(f"[Utils] Can't visualize. Basis must be a {2}x{2} matrix, but got a {basis.nrows()}x{basis.ncols()} one")
            thread = Utils.start_visualization(basis, closest_vector, point, title="Babai rounding technique")
            thread.join()
        
        return closest_vector