import subprocess, msprime, tskit, numpy, sys, random


if __name__ == "__main__":
    ts_path = sys.argv[1] + ".5indv.trees"
    print("Input tree sequence: %s\n" % ts_path)
    seed = int(sys.argv[2]) #Needed to sample the same individuals
    print("Random seed: %i\n" % seed)
    numpy.random.seed(seed)
    wsize = int(sys.argv[3])
    print("Window size: %i\n" % wsize)

    ts = tskit.load(ts_path)

    # Measure the tree height at each base position (adapted from the SLiM manual, section 17.4 Detecting the “dip in diversity”: analyzing tree heights in Python
    height_for_pos = numpy.zeros(int(ts.sequence_length))
    for tree in ts.trees():
        mean_height = numpy.mean([tree.time(root) for root in tree.roots])
        left, right = map(int, tree.interval)
        height_for_pos[left: right] = mean_height

    # Average per window and write results:
    outfile = sys.argv[1] + ("_w%i" % wsize) + ".csv"
    with open(outfile, "w") as output:
        output.write("Start,Stop,AverageTmrca\n")
        i = 0
        while i < int(ts.sequence_length):
            j = min(i + wsize, int(ts.sequence_length))
            mean_tmrca = numpy.mean(height_for_pos[i: j])
            output.write("%i,%i,%f\n" % (i, j, mean_tmrca))
            i = i + wsize

