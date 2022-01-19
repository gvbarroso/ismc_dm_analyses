import subprocess, msprime, pyslim, numpy, pandas, sys, random


if __name__ == "__main__":
    ts_path = sys.argv[1] + ".trees"
    print("Input tree sequence: %s\n" % ts_path)
    seed = int(sys.argv[2]) #Needed to sample the same individuals
    print("Random seed: %i\n" % seed)
    numpy.random.seed(seed)
    wsize = int(sys.argv[3])
    print("Window size: %i\n" % wsize)

    orig_ts = pyslim.load(ts_path)

    positions = []
    rates = []
    with open('../Comeron_tables/Comeron_100kb_chr2L.txt', 'r') as file:
        for line in file:
            components = line.split("\t")
            positions.append(float(components[0]))
            rates.append(1e-8 * float(components[1]))

    # step 1
    positions.insert(0, 0)
    # step 2
    rates.append(0.0)
    # step 3
    positions[-1] = 23.51e6
    positions[-1] += 1

    recomb_map = msprime.RecombinationMap(positions, rates)
    rts = orig_ts.recapitate(recombination_map = recomb_map, Ne = 10000)

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



    # Measure the tree height at each base position (adapted from the SLiM manual, section 17.4 Detecting the “dip in diversity”: analyzing tree heights in Python
    height_for_pos = numpy.zeros(int(sts.sequence_length))
    for tree in sts.trees():
        mean_height = numpy.mean([tree.time(root) for root in tree.roots])
        left, right = map(int, tree.interval)
        height_for_pos[left: right] = mean_height

    # Average per window and write results:
    outfile = sys.argv[1] + ("_w%i" % wsize) + ".csv"
    with open(outfile, "w") as output:
        output.write("Start,Stop,AverageTmrca\n")
        i = 0
        while i < int(sts.sequence_length):
            j = min(i + wsize, int(sts.sequence_length))
            mean_tmrca = numpy.mean(height_for_pos[i: j])
            output.write("%i,%i,%f\n" % (i, j, mean_tmrca))
            i = i + wsize

