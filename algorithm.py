from math import sqrt, floor, exp, log
import itertools

def legendre(a, p):
    return pow(a, (p - 1) // 2, p)

def tonelli_shanks_algo(n, p):
    """"Citation: This algorithm implementation was beyond scope for me to implement so I'm utilizing the implementation from https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Python"""

    assert legendre(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    z = 1
    for itz in range(2, p):
        z = itz
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        i = 0
        for i2 in range(1, m):
            i = i2
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r


def gaussian_elimination(matrix):
    rows, cols = len(matrix), len(matrix[0])
    rref_matrix= []

    for pivot_row in range(rows):
        # Find the pivot
        pivot = None
        for j in range(cols):
            if matrix[pivot_row][j] == 1:
                pivot = j
                break
        
        # If pivot is not found, skip this row
        if pivot is None:
            continue
        
        # Make the pivot row have 1 as the pivot
        for j in range(cols):
            matrix[pivot_row][j] = (matrix[pivot_row][j] + matrix[pivot_row][pivot]) % 2
        
        # Perform row operations to eliminate other entries in the pivot column
        for i in range(rows):
            if i != pivot_row and matrix[i][pivot] == 1:
                for j in range(cols):
                    matrix[i][j] = (matrix[i][j] + matrix[pivot_row][j]) % 2
    
    # Back substitution
    for i in range(rows - 1, -1, -1):
        pivot = None
        for j in range(cols):
            if matrix[i][j] == 1:
                pivot = j
                break
        if pivot is not None:
            row = [0] * cols
            row[pivot] = matrix[i][-1]
            for k in range(i):
                if matrix[k][pivot] == 1:
                    matrix[k][-1] = (matrix[k][-1] + matrix[k][pivot]) % 2
            rref_matrix.append(row)

    return rref_matrix

def isPrime(n):
    # Miller-Rabin primality test
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    even = n - 1
    k = 0
    q = even
    while q % 2 == 0:
        q //= 2
        k += 1
    test_amount = 25
    if test_amount > n - 1: 
        test_amount = n - 1
    for a in range(1, test_amount):
        if euclidean_algorithm(a, n) == 1:
            if (pow(a, q, n) == 1):
                continue
            value_to_check = -1 % n
            value_found = False
            for r in range(k):
                if (pow(a, (2**r) * q, n) == value_to_check):
                    value_found = True
                    break
            if value_found == False:
                return False
    return True

def insert_into_list(list, value):
    for i in range(len(list)):
        if value < list[i]:
            list.insert(i, value)
            return
    list.append(value)

def euclidean_algorithm(a, b):
    while b:
        a, b = b, a % b
    return abs(a)

def keep_only_indices(list, indices):
    new_list = []
    for i in indices:
        new_list.append(list[i])
    return new_list

def sieve_algorithm(N, B, X):
    #produces a list of as many of the a as possible and a list of a squared factored into small primes (so a list of lists)
    def F(T):
        return T**2 - N

    # else, X is a custom value
    list_of_a = []
    list_of_c = []
    list_of_c_factored = []

     # Find all of the primes up to B, can be done with sieve of Eratosthenes
    print("Finding primes up to B...")
    factor_base = [2]
    for i in range(3, B):
        is_prime = isPrime(i)
        if is_prime:
            if legendre(N, i) == 1:
                factor_base.append(i)

    
    a = floor(sqrt(N)) + 1 
    # I'm including 2 default options for X
    if X == 0:
        # This bound is from (Landquis, 2001), it is what he calls the sieving interval. It also changes the bottom bound of the interval
        b = sqrt(N) + exp(sqrt(log(N) * log(log(N))))**(3 * sqrt(2)/4)
        a = sqrt(N) - exp(sqrt(log(N) * log(log(N))))**(3 * sqrt(2)/4)
    elif X==1:
        b = 6 * len(factor_base) + a
    else:
        b = X * len(factor_base) + a
    

    # # Calculate F(curr) for all values in the range [a, b], naive approach
    # while (curr <= b):
    #     curr_val = F(curr)
    #     list_of_a.append(curr)
    #     list_of_c.append(curr_val)
    #     curr += 1
    
    # Calculate F(curr) for values that are divisible by primes in factor base using Tonelli-Shanks
    print("Calculating F(T) with Tonelli-Shanks...")
    for prime in factor_base:
        if prime == 2:
            alpha_p = 0
            beta_p = 0
            i = a
        else:      
            alpha_p = tonelli_shanks_algo(N, prime)
            beta_p = prime - alpha_p % prime
            i = floor((a - alpha_p) / prime)
            if floor(a - beta_p) > i:
                i = floor((a - beta_p) / prime)
        # i is the index of the first value above a
        while (alpha_p + i * prime < b) and (beta_p + i * prime < b):
            if alpha_p + i * prime not in list_of_a:
                list_of_a.append(alpha_p + i * prime)
                list_of_c.append(F(alpha_p + i * prime))
            if beta_p + i * prime not in list_of_a:
                list_of_a.append(beta_p + i * prime)
                list_of_c.append(F(beta_p + i * prime))
            i += 1   
    
    #primes_with_powers = factor_base.copy()
    # Introduce all of the prime powers, not used in the current implementation
    # for prime in factor_base:
    #     curr = prime
    #     while curr**2 < B:
    #         #insert_into_list(factor_base, curr)
    #         #appending does the same thing
    #         curr = curr**2
    #         primes_with_powers.append(prime)

    # Divide by all of the small prime factors, track the indices of values that are 1. Additionally, we keep track of the factors that divide the numbers
    print("Performing sieve...")
    indices_of_ones = []
    i = 0
    all_factor_occurences = []
    for elt in list_of_c:
        factors = []
        factor_occurences = []
        # factor_occurences is used for the linear algebra part, as it tracks the power of each prime in the factorization
        t = 0
        for prime in factor_base:
            factor_occurences.append(0)
            #running_prime_power = prime
            while elt % prime == 0:
                # if running_prime_power < B:
                # we do while to account for prime powers
                factor_occurences[t] = (factor_occurences[t] + 1) % 2
                factors.append(prime)
                elt = elt / prime
            t += 1
        if elt == 1:
            print("Product of small factors found!")
            indices_of_ones.append(i)
            list_of_c_factored.append(factors)
        i += 1
        all_factor_occurences.append(factor_occurences)
    
    # Now, we remove everything else from the original list of a values and the list of factors

    list_of_a = keep_only_indices(list_of_a, indices_of_ones)
    all_factor_occurences = keep_only_indices(all_factor_occurences, indices_of_ones)

    if list_of_a == []:
        print("No factors found. Try a higher B value or X value.")
        return None
    elif len(list_of_a) < len(factor_base):
        print("Not enough factors found. Try a higher B value or X value.")
        return None
    print(list_of_c_factored)
    return list_of_a, list_of_c_factored, factor_base, all_factor_occurences

def transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

def unique(list):
    unique_list = []
    for x in list:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

def generate_solutions(A, b):
    # Convert to row echelon form
    augmented_matrix = [row + [value] for row, value in zip(A, b)]
    rref_matrix = gaussian_elimination(augmented_matrix)

    # Find the free variables
    free_vars = []
    for i in range(len(rref_matrix[0]) - 1):
        column = [row[i] for row in rref_matrix]
        if sum(column) == 0:
            free_vars.append(i)

    row_count = len(rref_matrix)
    column_count = len(rref_matrix)
    # Find the constants
    constants = [row[-1] for row in A[:row_count]]  

    # Generate all possible solutions
    for free_vals in itertools.product([0, 1], repeat=len(free_vars)):
        solution = [0] * len(A[0])
        for i, value in enumerate(free_vals):
            solution[free_vars[i]] = value
        for i in range(len(constants)):
            for j in range(len(free_vars)):
                constants[i] = (constants[i] + A[i][free_vars[j]] * solution[free_vars[j]]) % 2
            solution[row_count - len(constants) + i] = constants[i]        
        yield solution

def quadratic_sieve_factorize(N, B = 0, X = 0):
    # Step 1: Relation building
    #How can we efficiently find many numbers a > √N such that each a^2 (mod N) is B-smooth?
    print("Building relations:")
    # I put this in as an option to allow for custom B values
    if B == 0:
        B = round(exp(sqrt(log(N) * log(log(N))))**(1/sqrt(2)))
    a_list, c_factored, primes, primes_vectors = sieve_algorithm(N, B, X)

    #step 2: Elimination
    # Now, we take the list of c factored, and we try the different products with those elements until we get every number appearing evenly
    # If each factor in the product of those cs appears an even amount of times, then the product is a perfect square, the square root of that is b and a is the product of all of the original as that formed the cs used to make b
    print("Performing elimination:")
    primes_vectors = transpose(primes_vectors)

    zero_array = [0] * len(primes_vectors)
    
    for u_vector in generate_solutions(primes_vectors, zero_array):
        u_vector = u_vector
        # Now we have the u vector, we can multiply the corresponding a values together to get b
        a = 1
        for i in range(len(u_vector)):
            if u_vector[i] == 1:
                a *= a_list[i]
        
        b = 1
        for i in range(len(u_vector)):
            if u_vector[i] == 1:
                for factor in c_factored[i]:
                    b *= factor
        #step 3: GCD computation
        #p = gcd(N, a − b), this is the factor
        # q = N/p

        print("Calculating GCD for candidate factors:")
        p = euclidean_algorithm(N, a - b)
        q = N // p
        assert(p * q == N)
        if p != 1 and q != 1:
            return p, q
        else:
            print("No factors found from this pair")
        
    print("No factors found :( Try a higher B or X value)")
    return None

print(quadratic_sieve_factorize(360379, 50, 150))
print("Expected Factors: 557, 647")
print(quadratic_sieve_factorize(914387, 100, 0))
print("Expected Factors: 829, 1103")
print(quadratic_sieve_factorize(123456789, 100))
print("Expected Factors: 3, 41152263")

#This is about the most it got in a reasonable amount of time, which is 47 bits
print(quadratic_sieve_factorize(123456789101112, 0, 0))
print("Expected Factors: 8, 15432098637639")





# This was my attempt to make a functioning more efficient elimination step, but I just couldn't get it to work
def generate_basis_solutions(matrix):
    rows, cols = len(matrix), len(matrix[0])
    rref_matrix= []
    # augment the matrix with the 0 vector as the last column
    for i in range(rows):
        matrix[i].append(0)
    rows = len(matrix)
    cols = len(matrix[0])

    # Forward elimination
    for col in range(min(rows, cols)):
        # Find pivot row
        pivot_row = col
        while pivot_row < rows and matrix[pivot_row][col] == 0:
            pivot_row += 1

        if pivot_row == rows:
            # No non-zero pivot in this column, move to the next column
            continue

        # Swap pivot row with current row (if needed)
        if pivot_row != col:
            matrix[col], matrix[pivot_row] = matrix[pivot_row], matrix[col]

        # Make pivot element 1
        pivot_value = matrix[col][col]
        for j in range(col, cols):
            matrix[col][j] /= pivot_value

        # Eliminate non-zero elements below the pivot
        for i in range(col + 1, rows):
            factor = matrix[i][col]
            for j in range(col, cols):
                matrix[i][j] -= factor * matrix[col][j]

    # Back substitution
    for col in range(min(rows, cols) - 1, 0, -1):
        for row in range(col - 1, -1, -1):
            factor = matrix[row][col]
            for j in range(col, cols):
                matrix[row][j] -= factor * matrix[col][j]

    rref_matrix = matrix
    
    # Now find each pivot again
    free_var_cols = []
    pivot_indices = []

    for row in matrix:
        pivot_index = next((i for i, x in enumerate(row) if x != 0), None)
        if pivot_index is not None:
            pivot_indices.append(pivot_index)

    for col in range(cols):
        if col not in pivot_indices and col != cols - 1:
            free_var_cols.append(col)
    
    for i in range(rows):
        for j in range(cols):
            rref_matrix[i][j] = rref_matrix[i][j] % 2

    #set all free variable columns to 0 except 1 and back substitute to find solution, then do for all free variable columns to get basis vectors

    basis_solutions = []
    for free_var_col in free_var_cols:
        temp_matrix = rref_matrix.copy()
        # free_var is 1, so we don't touch anything in the column
        for i in free_var_cols:
            if i != free_var_col:
                for j in range(rows):
                    # set all the other free_vars to 0
                    temp_matrix[j][i] = 0
        # Now we have a matrix with no free variables, we can back substitute and find the coefficients in the last column
        for i in range(rows - 1, -1, -1):
            pivot = None
            for j in range(cols):
                if temp_matrix[i][j] == 1:
                    pivot = j
                    break
            if pivot is not None:
                row = [0] * cols
                row[pivot] = temp_matrix[i][-1]
                for k in range(i):
                    if temp_matrix[k][pivot] == 1:
                        temp_matrix[k][-1] = (temp_matrix[k][-1] + temp_matrix[k][pivot]) % 2
        # figure out if variables are 0 or 1
        # this is the part that is inaccurate, along with the back substitution, I cannot figure out how to use the results of the back substitution to find the values of the basis vectors. This code at least manages to get the length of the basis I think
        basis_vector = []
        for j in range(cols):
            sum = 0
            for k in range(rows):
                sum += temp_matrix[k][j]
            if sum == 0:
                basis_vector[:-1].append(0)
            else:
                basis_vector.append(1)

        basis_solutions.append(basis_vector)

    return basis_solutions



# basis = generate_basis_solutions([[1,0,1,0,0,1,1,1,1,1,0,0,1,1,1,1,0,0,0,1], 
# [0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,1,1],
# [1,0,0,1,1,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0],
# [0,0,1,0,0,1,0,0,1,1,1,1,0,1,0,0,0,0,0,1],
# [1,1,0,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0],
# [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
# [0,1,0,0,0,1,1,1,0,0,0,1,0,1,1,0,1,1,0,0],
# [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0],
# [1,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0],
# [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],
# [0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0],
# [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1],
# [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
# [0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0],
# [0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0]])

# print(basis)
# print(len(basis))