# script for profiling halLiftover
import numpy as np
import time
from argparse import ArgumentParser
from subprocess import check_call

NUM_SAMPLES = 5

parser = ArgumentParser()
parser.add_argument('hal', default="/mnt/ntfs/RC_checkOutgroupInclusion_limit500k_10iter.hal")
parser.add_argument('--no-perf', action='store_true')
parser.add_argument('--no-collapse', action='store_true', help="don't collapse self-recursive functions")
opts = parser.parse_args()

times = []
for i in range(NUM_SAMPLES):
    check_call("echo 3 > /proc/sys/vm/drop_caches", shell=True)
    startTime = time.time()
    check_call(["halLiftover", opts.hal, "human", "human.bed", "mouse", "mouse.bed"])
    stopTime = time.time()
    times.append(stopTime-startTime)

times = np.array(times)
print np.mean(times), "+/-", np.std(times, ddof=1)

if not opts.no_perf:
    # Another run to generate a flamegraph. This is separate from the
    # other runs to avoid any potential effects of perf overhead
    # (though they should be minuscule).
    check_call("perf record -F 612 -g --call-graph dwarf -- halLiftover %s human human.bed mouse mouse.bed" % opts.hal, shell=True)
    check_call("perf script | stackcollapse-perf.pl > halLiftover.stacks", shell=True)
    with open('halLiftover.stacks') as input, open('halLiftover.stacks.merged', 'w') as output:
        for line in input:
            stacks = line.strip().split(';')
            try:
                i = stacks.index('hal::BlockLiftover::liftInterval')
            except ValueError:
                pass
            else:
                stacks = stacks[i:]
            if not opts.no_collapse:
                # Collapse recursive runs of the same function into a single instance.
                new_stacks = []
                prev_entry = None
                for entry in stacks:
                    if entry != prev_entry:
                        new_stacks.append(entry)
                    prev_entry = entry
                stacks = new_stacks
            output.write(";".join(stacks) + "\n")
    check_call("flamegraph.pl halLiftover.stacks.merged > halLiftover.svg", shell=True)
