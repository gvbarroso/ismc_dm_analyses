import subprocess, msprime, tskit, numpy, pandas, sys, random, warnings


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
    ts_path = sys.argv[1] + ".5indv.trees"
    print("Input tree sequence: %s\n" % ts_path)
    seed = int(sys.argv[2])
    print("Random seed: %i\n" % seed)
    numpy.random.seed(seed)

    ts = tskit.load(ts_path)

    mutation_model = msprime.SLiMMutationModel(0, 0)

    mts = msprime.sim_mutations(ts, rate = 1e-7, model = mutation_model, keep = True, random_seed = seed)

    t = mts.tables
    t.sites.clear()
    t.mutations.clear()


    for s in mts.sites():
        alleles = make_alleles(len(s.mutations) + 1)
        t.sites.append(s.replace(ancestral_state=alleles[0]))
        for j, m in enumerate(s.mutations):
            t.mutations.append(m.replace(derived_state=alleles[j+1][0])) # If there are more than 4 states, takes only first character fo 5th state 

    mts = t.tree_sequence()

    outfile = "Homogeneous/" + sys.argv[1] + "_5.vcf"

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

    with open("Homogeneous/" + sys.argv[1] + "_5.fasta", 'w') as output:
        for i, s in enumerate(S) :
            output.write(">hap%i\n" % i)
            output.write(s.tobytes().decode(encoding))
            output.write("\n")


