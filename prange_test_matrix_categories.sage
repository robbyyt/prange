MAX_ITERATIONS = 100000

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
        
        if iteration_count > MAX_ITERATIONS:
            raise Exception("Maximum iterations exceeded in inner loop")

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
            if verbose:
                print("ISD finished after %d iterations of outer loop and an average of %d iterations of inner loop"
                      % (outer_iter_count, mean(inner_iter_counts)))
            
            return matrix(e_curr) * P.T, outer_iter_count, mean(inner_iter_counts)
        
        if outer_iter_count > MAX_ITERATIONS:
            raise Exception("Maximum iterations exceeded in outer loop")


def randvect(size=4):
    return vector([1 if random() >= 0.5 else 0 for t in range(size)])

def test_prange(H, order):
    for iter_alg in range(1):
        print("Iteration %d\n" % iter_alg)
        inner_it_avgs = []
        outer_its = []
        for curr_weight in range(1, order):
            print("Runnig alg for syndrome weight %d" % curr_weight)
            for iter_per_weight in range(1, 3):
                syndrome = matrix(order, 1, randvect(order))
                try:
                    e, outer_it_count, inner_it_avg = prange_ISD(H, curr_weight, syndrome)
                except Exception as err:
                    print("Failed to find a solution!\n")
                    print("Syndrome:\n", syndrome)
                    print("Weight:\n", curr_weight)
                    print("---------")
                    print(err)
                    continue

                if H * e.transpose() == syndrome:
                    print("ISD returned correct result")
                    inner_it_avgs.append(inner_it_avg)
                    outer_its.append(outer_it_count)
            
            print("ISD ran for 2 iterations on weight %d with an average of %d outer loop runs and an average of %d iterations of inner loop"
                      % (curr_weight, mean(outer_its), mean(inner_it_avgs)))
            print("--------------------------------------------------------------")



ORDER = 10
C = codes.HammingCode(GF(2), ORDER)
PC = C.parity_check_matrix()

test_prange(PC, ORDER)