import subprocess, msprime, pyslim, numpy, pandas, sys, random


# see: https://github.com/tskit-dev/pyslim/issues/168
def make_alleles(n):
    nucs = ['A', 'C', 'G', 'T']
    alleles = nucs.copy()
    numpy.random.shuffle(alleles)
    while len(alleles) < n:
        new = [a + b for a in alleles for b in nucs]
        numpy.random.shuffle(new)
        alleles.extend(new)
    return alleles[:n]

if __name__ == "__main__":
    ts_path = sys.argv[1] + ".trees"
    print("Input tree sequence: %s\n" % ts_path)
    seed = int(sys.argv[2])
    print("Random seed: %i\n" % seed)
    numpy.random.seed(seed)

    orig_ts = pyslim.load(ts_path)

    rec_positions = []
    rec_rates = []
    with open('../Comeron_tables/Comeron_100kb_chr2L.txt', 'r') as file:
        for line in file:
            components = line.split("\t")
            rec_positions.append(float(components[0]))
            rec_rates.append(1e-7 * float(components[1])) #10r

    # step 1
    rec_positions.insert(0, 0)
    # step 2
    rec_rates.append(0.0)
    # step 3
    rec_positions[-1] = 23.51e6
    rec_positions[-1] += 1

    recomb_map = msprime.RecombinationMap(rec_positions, rec_rates)
    rts = recapitate(orig_ts, recombination_map = recomb_map, ancestral_Ne = 10000)

    assert(max([t.num_roots for t in rts.trees()]) == 1)


    alive_inds = rts.individuals_alive_at(0)
    keep_indivs = numpy.random.choice(alive_inds, 5, replace = False)
    print(f"Individuals kept: {keep_indivs}")
    
    keep_nodes = []
    for i in keep_indivs:
        keep_nodes.extend(rts.individual(i).nodes)
    sts = rts.simplify(keep_nodes)
    print(f"Nodes kept: {keep_nodes}")

    print(f"Before, there were {rts.num_samples} sample nodes (and {rts.num_individuals} individuals) "
          f"in the tree sequence, and now there are {sts.num_samples} sample nodes "
          f"(and {sts.num_individuals} individuals).")

    mutation_model = msprime.SLiMMutationModel(0, 0)

    mut_positions = []
    mut_rates = []
    with open('../MutationMap.csv', 'r') as file:
        next(file) # Skip header
        for line in file:
            components = line.split(",")
            mut_positions.append(int(components[0]))
            mut_rates.append(1e-7 * float(components[1]))
    mut_positions.append(23.51e6)
    mut_positions[-1] += 1

    mutation_map = msprime.RateMap(position = mut_positions, rate = mut_rates)

    mts = msprime.sim_mutations(sts, rate = mutation_map, model = mutation_model, keep = True, random_seed = seed)

    t = mts.tables
    t.sites.clear()
    t.mutations.clear()


    for s in mts.sites():
        alleles = make_alleles(len(s.mutations) + 1)
        t.sites.append(s.replace(ancestral_state=alleles[0]))
        for j, m in enumerate(s.mutations):
            t.mutations.append(m.replace(derived_state=alleles[j+1][0])) # If there are more than 4 states, takes only first character fo 5th state 

    mts = t.tree_sequence()

    outfile = "NonHomogeneous/" + sys.argv[1] + "_5.vcf"

    with open(outfile, "w") as vcf_file:
        mts.write_vcf(vcf_file) #, position_transform = "legacy") # 'legacy' is used to avoid SNPs at the same position

    # Now print sequences as fasta (adapted from https://github.com/tskit-dev/tskit/issues/353#issuecomment-745275431)
    # Non-variable sites are all A with mutations as T.

    encoding = "ascii"
    S = numpy.zeros((mts.num_samples, int(mts.sequence_length)), dtype = numpy.int8)
    pad_value = ord("A".encode(encoding))
    S[:] = pad_value

    for variant in mts.variants():
        alleles = numpy.array([ord(allele.encode(encoding)) for allele in variant.alleles], dtype = numpy.int8)
        coord = int(variant.site.position)
        if S[0, coord] == pad_value:
            # This sites hasn't been seen before, so initialise it to the
            # ancestral state
            S[:, coord] = alleles[0]
        # Update the site with any derived states
        derived_samples = variant.genotypes > 0
        S[derived_samples, coord] = alleles[variant.genotypes[derived_samples]]

    # print out the sequences

    with open("NonHomogeneous/" + sys.argv[1] + "_5.fasta", 'w') as output:
        for i, s in enumerate(S) :
            output.write(">hap%i\n" % i)
            output.write(s.tobytes().decode(encoding))
            output.write("\n")


