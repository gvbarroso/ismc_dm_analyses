import subprocess, msprime, pyslim, numpy, sys, warnings


if __name__ == "__main__":
    ts_path = sys.argv[1] + ".trees"
    print("Input tree sequence: %s\n" % ts_path)
    seed = int(sys.argv[2])
    print("Random seed: %i\n" % seed)
    numpy.random.seed(seed)

    orig_ts = pyslim.load(ts_path)

    positions = [int(x * 1e6) for x in [0, 5, 10, 15, 20, 25, 30]]
    positions[-1] = positions[-1] + 1
    rates = [x * 1e-7 for x in [0.1, 1.5, 0.1, 1.5, 0.1, 1.5]] #10r
    recomb_map = msprime.RateMap(position = positions, rate = rates)
    warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning) #If tree sequence is in older format
    rts = pyslim.recapitate(orig_ts, recombination_rate = recomb_map, ancestral_Ne = 10000)

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


    sts.dump(sys.argv[1] + ".5indv.trees")

