SAMPLE_PC_MATRIX = matrix([
    [1,0,1,0,1,0,1],
    [0,1,1,0,0,1,1],
    [0,0,0,1,1,1,1]
    ])


SAMPLE_SYNDROME = matrix(3,1,[1,0,1])

SAMPLE_WEIGHT= 2

SAMPLE_PC_MATRIX_2 = Matrix(ZZ, 
    matrix([
    [0,0,1,0,1,0,0,0,0,0,0,1],
    [0,1,1,1,0,0,0,0,1,0,0,0],
    [0,0,0,0,1,1,0,1,0,1,1,1],
    [1,1,0,0,0,1,1,1,0,0,0,0],
    [1,1,1,1,0,0,1,0,1,1,1,1],
    [0,0,1,0,0,0,0,0,1,1,1,0],
    [0,1,0,0,0,0,0,1,0,1,0,1],
    [0,1,1,1,1,0,1,1,1,0,0,1]
    ])
)

SAMPLE_SYNDROME_2 = Matrix(ZZ, matrix(8,1,[0,1,1,0,0,1,0,0]))
SAMPLE_WEIGHT_2 = 3

def generateRandomPermutationMatrix(n):
    return Permutations(n).random_element().to_matrix()

def is_of_desired_weight(arr, weight):
    w_curr = 0
    for x in arr:
        if x != 0:
            w_curr += 1
    
    if w_curr != weight:
        return False

    return True

def prange_inner_loop(H):
    r, n = H.dimensions()
    # inner repeat-until block
    iteration_count = 0
    while True:
        P = generateRandomPermutationMatrix(n)
        H_curr = H * P
        H_curr = block_matrix(1, (H_curr, identity_matrix(r)))
        RRE = H_curr.rref()
        U = RRE[:, n:n+r] 
        V = RRE[:, r:n]
        W = RRE[:, 0:r]
        # print("Inner loop Iteration: %d\n" % (iteration_count + 1))
        # print("W:\n", W, '\n-----')
        # print("U:\n", U, '\n-----')
        # print("V:\n", V, '\n-----')
        iteration_count += 1
        if W == identity_matrix(r):
            return P,V,U, iteration_count

def prange_ISD(H, t, s, verbose=False):
    inner_iter_counts = []
    outer_iter_count = 0
    while True:
        P, V, U, inner_iteration_count = prange_inner_loop(H)
        inner_iter_counts.append(inner_iteration_count)
        
        r,n = H.dimensions()
        s_curr = U * s
        s_curr = s_curr.transpose().list()
        e_curr =  s_curr + [0] * (n - len(s_curr))
        
        outer_iter_count += 1
        
        if is_of_desired_weight(e_curr, t):
            print("e: ", e_curr)
            print(P.transpose())
            if verbose:
                print("ISD finished after %d iterations of outer loop and an average of %d iterations of inner loop"
                      % (outer_iter_count, mean(inner_iter_counts)))
            
            return matrix(e_curr) * P.T


e = prange_ISD(SAMPLE_PC_MATRIX_2, SAMPLE_WEIGHT_2, SAMPLE_SYNDROME_2, verbose=True)
print(e)
print("OUTPUT:\n", SAMPLE_PC_MATRIX_2 * e.transpose())
