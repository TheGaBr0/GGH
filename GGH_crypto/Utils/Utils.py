from flint import fmpz_mat, fmpq_mat, fmpq, fmpz
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal, getcontext
import math
import os
import subprocess
import ast
from fractions import Fraction
import re
import time
import numpy as np
import random

class Utils:
    """
    A utility class providing various helper methods for lattice-based cryptography operations.

    This class includes methods for linear algebra operations, lattice visualizations,
    mathematical computations, and file I/O operations specific to lattice-based cryptography.

    Class Methods:
        gram_schmidt(basis): Performs Gram-Schmidt orthogonalization on a given basis.
        visualize_lattice(basis_1, basis_1_cvp, point, basis_2=None, basis_2_cvp=None, title="Lattice Plot", limit=5):
            Visualizes a 2D lattice with given bases and points.
        embedding(basis, ciphertext, visualize=False, GGH=False, BKZ=False, block=20, pruned=False, precision=100, bkzautoabort=True, bkzmaxloops=None, nolll=False):
            Performs lattice embedding for closest vector problem (CVP) solving.
        npsp_to_fmpq_mat(basis): Converts a numpy array or sympy matrix to an fmpq_mat.
        npsp_to_fmpz_mat(basis): Converts a numpy array or sympy matrix to an fmpz_mat.
        vector_l1_norm(row): Calculates the L1 norm of a vector.
        vector_l2_norm(row): Calculates the L2 norm of a vector.
        get_hadamard_ratio(basis=None, precision=10): Calculates the Hadamard ratio of a basis.
        write_matrix_to_file(matrix, filename): Writes a matrix to a file.
        load_matrix_from_file(filename, matrix_type='fmpq'): Loads a matrix from a file.
        BKZ_reduction(matrix, block=20, pruned=False, precision=90, bkzautoabort=True, bkzmaxloops=None, nolll=False):
            Performs BKZ (Block Korkine-Zolotarev) lattice reduction.
        babai_rounding(basis, point, visualize=False): Performs Babai's rounding algorithm for CVP.

    Note:
        This class assumes the availability of certain libraries like flint, matplotlib, and numpy.
        It also interacts with system commands, particularly for BKZ reduction using fplll.
    """
    def gram_schmidt(basis):
        """
        Performs Gram-Schmidt orthogonalization on a given basis.

        Args:
            basis (numpy.ndarray): The input basis.

        Returns:
            numpy.ndarray: The orthogonalized basis.
        """
        ortho_basis = basis[0:1,:].copy()
        for i in range(1, basis.shape[0]):
            proj = np.diag((basis[i,:].dot(ortho_basis.T)/np.linalg.norm(ortho_basis,axis=1)**2).flat).dot(ortho_basis)
            ortho_basis = np.vstack((ortho_basis, basis[i,:] - proj.sum(0)))
        return ortho_basis
    
    def visualize_lattice(basis_1, basis_1_cvp, point, basis_2=None, basis_2_cvp=None, title="Lattice Plot", limit=5):
        """
        Visualizes a 2D lattice with given bases and points.

        Args:
            basis_1 (fmpz_mat): The first basis.
            basis_1_cvp (fmpz_mat): The closest vector point for basis_1.
            point (fmpz_mat): The target point.
            basis_2 (fmpz_mat, optional): The second basis.
            basis_2_cvp (fmpz_mat, optional): The closest vector point for basis_2.
            title (str, optional): The title of the plot.
            limit (int, optional): The limit for the lattice points to plot.
        """
        basis_1_np = np.array(basis_1.tolist()).astype(int)

        plt.rcParams.update({'font.size': 14})

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
        
        fig, ax = plt.subplots(figsize=(10, 10))
        
        ax.plot(lattice_points[:, 0], lattice_points[:, 1], 'bo', color='black')
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_title(title)
        
        # Add Cartesian axes
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.axvline(x=0, color='k', linewidth=0.5)
        
        label_w = "w"
        label_w_2 = "w_2"
        label_t = "t"
        
        ax.scatter(point_np[0][0], point_np[0][1])
        
        if basis_2_cvp is not None:
            ax.scatter(int(basis_2_cvp[0, 0]), int(basis_2_cvp[0, 1]), color='red', s=70)
            
        if basis_1_cvp is not None:
            ax.scatter(int(basis_1_cvp[0, 0]), int(basis_1_cvp[0, 1]), color='blueviolet', s=70)
            
        ax.annotate(label_t, 
                    (point_np[0][0], point_np[0][1]), 
                    textcoords="offset points", 
                    xytext=(10, -5), 
                    ha='center')
                    
        if basis_1_cvp is not None:
            ax.annotate(label_w, 
                        (int(basis_1_cvp[0, 0]), int(basis_1_cvp[0, 1])), 
                        textcoords="offset points",
                        xytext=(0, 10), 
                        ha='center')
                        
        if basis_2_cvp is not None:
            ax.annotate(label_w_2, 
                        (int(basis_2_cvp[0, 0]), int(basis_2_cvp[0, 1])), 
                        textcoords="offset points", 
                        xytext=(0, -12), 
                        ha='center') 
            
        ax.arrow(0, 0, basis_1_np[0][0], basis_1_np[0][1], lw=2, head_width=0.1, head_length=0.1, fc='blueviolet', ec='blueviolet')
        ax.arrow(0, 0, basis_1_np[1][0], basis_1_np[1][1], lw=2, head_width=0.1, head_length=0.1, fc='blueviolet', ec='blueviolet')
        
        if basis_2 is not None:
            ax.arrow(0, 0, basis_2_np[0][0], basis_2_np[0][1], lw=2, head_width=0.1, head_length=0.1, fc='red', ec='red')
            ax.arrow(0, 0, basis_2_np[1][0], basis_2_np[1][1], lw=2, head_width=0.1, head_length=0.1, fc='red', ec='red')  

        # Set equal aspect ratio
        ax.set_aspect('equal', adjustable='box')

        # Calculate the limits based on lattice points and basis vectors
        all_points = [lattice_points]
        
        # Ensure all arrays are 2D before appending
        all_points.append(basis_1_np.reshape(-1, 2))
        all_points.append(point_np.reshape(-1, 2))
        
        if basis_2 is not None:
            all_points.append(basis_2_np.reshape(-1, 2))
        if basis_1_cvp is not None:
            all_points.append(np.array(basis_1_cvp.tolist()).astype(int).reshape(-1, 2))
        if basis_2_cvp is not None:
            all_points.append(np.array(basis_2_cvp.tolist()).astype(int).reshape(-1, 2))

        all_points = np.vstack(all_points)

        x_min, y_min = np.min(all_points, axis=0)
        x_max, y_max = np.max(all_points, axis=0)

        # Add some padding
        padding = 0.5
        ax.set_xlim(x_min - padding, x_max + padding)
        ax.set_ylim(y_min - padding, y_max + padding)

        plt.tight_layout()
        plt.show()

    def embedding(basis, ciphertext, visualize=False, GGH=False, BKZ=False, block=20, 
                  pruned=False, precision=100, bkzautoabort=True, bkzmaxloops=None, nolll=False):
        """
        Performs lattice embedding for closest vector problem (CVP) solving.

        Args:
            basis (fmpz_mat): The lattice basis.
            ciphertext (fmpz_mat): The ciphertext vector.
            visualize (bool, optional): Whether to visualize the result.
            GGH (bool, optional): Whether to use GGH-specific constraints.
            BKZ (bool, optional): Whether to use BKZ reduction.
            block (int, optional): The block size for BKZ reduction.
            pruned (bool, optional): Whether to use pruning in BKZ.
            precision (int, optional): The precision for BKZ calculations.
            bkzautoabort (bool, optional): Whether to use auto-abort in BKZ.
            bkzmaxloops (int, optional): The maximum number of loops for BKZ.
            nolll (bool, optional): Whether to skip LLL in BKZ.

        Returns:
            fmpz_mat: The closest vector to the ciphertext in the lattice.
        """
        if basis.ncols() != ciphertext.ncols():
                raise ValueError(f"[Utils] Point is a {1}x{ciphertext.ncols()} matrix, but basis is a {basis.nrows()}x{basis.ncols()} one")
        
        n = basis.nrows()

        # Create the initial matrix with R and a column of zeros
        matrix_emb = fmpz_mat([[int(basis[i,j]) if j < n else 0 for j in range(n+1)] for i in range(n)])

        # Add t as the last row, with 1 as the last element
        last_row = [int(ciphertext[0,i]) for i in range(n)] + [1]
        matrix_emb = fmpz_mat(matrix_emb.tolist() + [last_row])
        
        if BKZ:
            matrix_emb = Utils.BKZ_reduction(matrix_emb, block=block, pruned=pruned, precision=precision,
                                   bkzautoabort=bkzautoabort, bkzmaxloops=bkzmaxloops, nolll=nolll)
        else:
            matrix_emb = matrix_emb.lll()

        min_norm = float('inf')
        shortest_vector = None

        for i in range(n + 1):
            i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n + 1)]])

            if GGH:
                if matrix_emb[i, n] == 1:  # Se l'ultimo elemento è 1
                    norm = Utils.vector_l2_norm(i_th_row)
                    if norm < min_norm:
                        min_norm = norm
                        # Prendi i primi n valori della riga più corta
                        shortest_vector = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])
            else:
                # Se GGH non è vero, considera comunque la riga per trovare il vettore più corto
                norm = Utils.vector_l2_norm(i_th_row)
                if norm < min_norm:
                    min_norm = norm
                    # Prendi i primi n valori della riga più corta
                    shortest_vector = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])
    
        # Se GGH è vero e nessun vettore con 1 come ultimo elemento è stato trovato,
        # ritorna comunque il vettore più corto trovato
        if GGH and shortest_vector is None:
            for i in range(n + 1):
                i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n + 1)]])
                norm = Utils.vector_l2_norm(i_th_row)
                if norm < min_norm:
                    min_norm = norm
                    # Prendi i primi n valori della riga più corta
                    shortest_vector = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])

        closest_vector = ciphertext - shortest_vector

        if visualize:
            if basis.nrows() != 2:
                raise ValueError(f"[Utils] Can't visualize. Basis must be a {2}x{2} matrix, but got a {basis.nrows()}x{basis.ncols()} one")
            Utils.visualize_lattice(basis, closest_vector, ciphertext, title="Embedding method")
            

        return closest_vector
    
    def npsp_to_fmpq_mat(basis):
        """
        Converts a numpy array or sympy matrix to an fmpq_mat.

        Args:
            basis (numpy.ndarray or sympy.Matrix): The input matrix.

        Returns:
            fmpq_mat: The converted matrix.
        """
        fractions = [[Fraction(item) for item in row] for row in basis.tolist()]
        return fmpq_mat([[fmpq(f.numerator, f.denominator) for f in row] for row in fractions])

    def npsp_to_fmpz_mat(basis):
        """
        Converts a numpy array or sympy matrix to an fmpz_mat.

        Args:
            basis (numpy.ndarray or sympy.Matrix): The input matrix.

        Returns:
            fmpz_mat: The converted matrix.
        """
        return fmpz_mat([[int(item) for item in sublist] for sublist in basis.tolist()])
        
    def vector_l1_norm(row):
        """
        Calculates the L1 norm of a vector.

        Args:
            row (fmpz_mat or fmpq_mat): The input vector.

        Returns:
            Decimal: The L1 norm of the vector.
        """
        if isinstance(row, fmpz_mat):
            row = fmpq_mat(row)
        getcontext().prec = 50
        return max(sum(abs(Decimal(int(x.numer())) / Decimal(int(x.denom()))) for x in row) for row in row.tolist())
        
    def vector_l2_norm(row):
        """
        Calculates the L2 norm of a vector.

        Args:
            row (fmpz_mat or fmpq_mat): The input vector.

        Returns:
            Decimal: The L2 norm of the vector.
        """
        if isinstance(row, fmpz_mat):
            row = fmpq_mat(row)
        getcontext().prec = 50
        return Decimal(sum((Decimal(int(x.numer())) / Decimal(int(x.denom()))) ** 2 for x in row)).sqrt()

    
    def get_hadamard_ratio(basis=None, precision=10):
        """
        Calculates the Hadamard ratio of a basis.

        Args:
            basis (fmpz_mat or fmpq_mat, optional): The input basis.
            precision (int, optional): The precision for calculations.

        Returns:
            tuple: (Decimal ratio, str formatted ratio)
        """
        norms = []
        dimension = basis.nrows()
        
        # Set a high precision for Decimal calculations
        getcontext().prec = precision
        
        for i in range(basis.nrows()):
            row = fmpz_mat([[basis[i, j] for j in range(basis.ncols())]])
            norm = Utils.vector_l2_norm(row)
            norms.append(Decimal(str(norm)))
        
        # Use log sum instead of direct multiplication
        log_denominator = sum(norm.ln() for norm in norms)
        log_numerator = abs(Decimal(basis.det().str())).ln()
        
        # Calculate the ratio using logs
        log_result = (log_numerator - log_denominator) / Decimal(dimension)
        result = log_result.exp()
        
        return result, f"{result:.{precision}f}"
    
    def write_matrix_to_file(matrix, filename):
        """
        Writes a matrix to a file.

        Args:
            matrix (fmpz_mat or fmpq_mat): The matrix to write.
            filename (str): The name of the file to write to.

        Returns:
            str: The full path of the written file.
        """
        filename = os.path.join(os.getcwd(), filename)

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
    
    def load_matrix_from_file(filename, matrix_type='fmpq'):
        """
        Loads a matrix from a file.

        Args:
            filename (str): The name of the file to read from.
            matrix_type (str, optional): The type of matrix to load ('fmpz' or 'fmpq').

        Returns:
            fmpz_mat or fmpq_mat: The loaded matrix.

        Raises:
            ValueError: If an invalid matrix_type is provided.
        """
        with open(os.path.join(os.getcwd(), filename), 'r') as file:
            content = file.read().replace(']', '],').replace(' ]', ']').replace(' ', ', ')[:-4] + ']'
    
        if matrix_type == 'fmpz':
            return fmpz_mat(ast.literal_eval(content))
        elif matrix_type == 'fmpq':
            data = [[Fraction(*map(int, elem.split('/') if '/' in elem else (elem, 1)))
                    for elem in re.findall(r'[-]?\d+(?:/\d+)?', row)]
                    for row in re.findall(r'\[(.*?)\]', content)]
            
            matrix = fmpq_mat(len(data), len(data[0]) if data else 0)

            for i, row in enumerate(data):
                for j, frac in enumerate(row):
                    matrix[i, j] = fmpq(frac.numerator, frac.denominator)
            return matrix
        else:
            raise ValueError("Invalid matrix_type. Use 'fmpz' or 'fmpq'.")

    def BKZ_reduction(matrix, block=20, pruned=False, precision=90, bkzautoabort=True, bkzmaxloops=None, nolll=False):
        """
        Performs BKZ (Block Korkine-Zolotarev) lattice reduction.

        Args:
            matrix (fmpz_mat): The input matrix to reduce.
            block (int, optional): The block size for BKZ.
            pruned (bool, optional): Whether to use pruning.
            precision (int, optional): The precision for calculations.
            bkzautoabort (bool, optional): Whether to use auto-abort.
            bkzmaxloops (int, optional): The maximum number of loops.
            nolll (bool, optional): Whether to skip LLL.

        Returns:
            tuple: (fmpz_mat reduced matrix, str error message or None)
        """
        input_path = Utils.write_matrix_to_file(matrix, f'input.txt')
        output_path = Utils.write_matrix_to_file(matrix, f'output.txt')

        if os.name == 'nt':
            command = f"wsl fplll input.txt -a bkz -b {block} -p {precision} -m wrapper -f mpfr"
            if pruned:
                command += " -s default.json"
            if bkzautoabort:
                command += " -bkzautoabort"
            if bkzmaxloops != None:
                command += f" -bkzmaxloops {bkzmaxloops}"
            if nolll:
                command += " -nolll"
            command += f" > output.txt"
        else:
            command = f"fplll input.txt -a bkz -b {block} -p {precision} -m wrapper -f mpfr"
            if pruned:
                command += " -s default.json"
            if bkzautoabort:
                command += " -bkzautoabort"
            if bkzmaxloops != None:
                command += f" -bkzmaxloops {bkzmaxloops}"
            if nolll:
                command += " -nolll"
            command += f" > output.txt"
        try:
            # Run the command and capture its output
            print(f"Reduction started with the following parameters:\n"
            f"  block: {block}\n"
            f"  pruned: {pruned}\n"
            f"  precision: {precision}\n"
            f"  bkzautoabort: {bkzautoabort}\n"
            f"  bkzmaxloops: {bkzmaxloops}\n"
            f"  nolll: {nolll}")
            print("Final command:\n"
                  f"{command}")
           
            time_now = time.time()
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, cwd=os.getcwd())
            output, error = process.communicate()
            output = output.decode('utf-8')
            error = error.decode('utf-8')

            if error:
                if str(error).strip() != "Failure: loops limit exceeded in BKZ":
                    print("Error during reduction:", error)
            else:
                error = None

            print(f"Reduction completed, time taken: {time.time() - time_now}")

            # Load the reduced matrix
            reduced_matrix = Utils.load_matrix_from_file(f"output.txt", "fmpz")

            os.remove(output_path)
            os.remove(input_path)

            return reduced_matrix, error

        except Exception as e:
            return None, str(e)
        
    def babai_rounding(basis, point, visualize=False):
        """
        Performs Babai's rounding algorithm for CVP.

        Args:
            basis (fmpz_mat): The lattice basis.
            point (fmpz_mat): The target point.
            visualize (bool, optional): Whether to visualize the result.

        Returns:
            fmpz_mat: The closest vector to the point in the lattice.

        Raises:
            ValueError: If visualization is requested for a non-2D lattice.
        """
        x = point * basis.inv()

        for i in range(x.nrows()):
            for j in range(x.ncols()):
                x[i,j] = round(x[i,j])
        
        closest_vector = x * basis 

        if visualize:
            if basis.nrows() != 2:
                raise ValueError(f"[Utils] Can't visualize. Basis must be a {2}x{2} matrix, but got a {basis.nrows()}x{basis.ncols()} one")
            Utils.visualize_lattice(basis, closest_vector, point, title="Babai rounding technique")
            
        return closest_vector
