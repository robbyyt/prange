MAX_ITERATIONS = 10_000
f = open("C:\\Users\\Robbz\\Documents\\MASTER-IMPLEMENTATIONS\\prange\\results\\hamming_order5", "a+")

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
        iteration_count += 1
        if W == identity_matrix(r):
            return P, V, U, iteration_count

        if iteration_count > MAX_ITERATIONS:
            raise Exception("Maximum iterations exceeded in inner loop")


def prange_ISD(H, t, s, verbose=False):
    inner_iter_counts = []
    outer_iter_count = 0
    while True:
        try:
            P, V, U, inner_iteration_count = prange_inner_loop(H)
        except:
            raise Exception("Exiting ISD algorithm, maximum iterations exceeded in inner loop")

        inner_iter_counts.append(inner_iteration_count)

        r, n = H.dimensions()
        s_curr = U * s
        s_curr = s_curr.transpose().list()
        e_curr = s_curr + [0] * (n - len(s_curr))

        outer_iter_count += 1

        if is_of_desired_weight(e_curr, t):
            if verbose:
                f.write("ISD finished after %d iterations of outer loop and an average of %d iterations of inner loop\n"
                      % (outer_iter_count, mean(inner_iter_counts)))

            return matrix(e_curr) * P.T, outer_iter_count, mean(inner_iter_counts)

        if outer_iter_count > MAX_ITERATIONS:
            raise Exception("Maximum iterations exceeded in outer loop")


def randvect(size=4):
    return vector([1 if random() >= 0.5 else 0 for t in range(size)])


def test_prange(H, order):
    for iter_alg in range(10):
        f.write("Iteration %d\n" % iter_alg)
        inner_it_avgs = []
        outer_its = []
        for curr_weight in range(1, order):
            f.write("Runnig alg for syndrome weight %d\n" % curr_weight)
            for iter_per_weight in range(10):
                syndrome = matrix(order, 1, randvect(order))
                try:
                    e, outer_it_count, inner_it_avg = prange_ISD(
                        H, curr_weight, syndrome)
                except Exception as err:

                    f.write("Failed to find a solution!\n")
                    f.write("Weight:\n %d" % curr_weight)
                    f.write("---------\n")
                    continue

                if H * e.transpose() == syndrome:
                    f.write("ISD returned correct result\n")
                    inner_it_avgs.append(inner_it_avg)
                    outer_its.append(outer_it_count)
            try:
                f.write("ISD ran for 2 iterations on weight %d with an average of %d outer loop runs and an average of %d iterations of inner loop\n"
                  % (curr_weight, mean(outer_its), mean(inner_it_avgs)))
                f.write("--------------------------------------------------------------\n")
            except:
                f.write("Could not print statistics!\n")


ORDER = 5
C = codes.HammingCode(GF(2), ORDER)
PC = C.parity_check_matrix()

# test_prange(PC, ORDER)

# F = GF(2 ^ 6)

# R.<x> = F[]

# g = x ^ 9  + 1

# L = [a for a in F.list() if g(a) != 0]

# C = codes.GoppaCode(g, L)
# PC = Matrix(QQ, C.parity_check_matrix())
# f.write(" DIMENSIONS: ", PC.dimensions())
# ORDER = PC.dimensions()[0]


# C = codes.GolayCode(GF(2), extended = False)
# PC = C.parity_check_matrix()
# ORDER = PC.dimensions()[0]

test_prange(PC, ORDER)
f.close()