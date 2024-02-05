#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import subprocess, argparse

def get_gflops(ts, mode, n):
    """Run matrix_prod program and report tile size and GFlop/s rate"""
    # Constuct command line and run command - GFlops rate is 4th token
    # in space-delimited output
    cmd = f'./matrix_prod -r {mode} {n} {n} {n} {ts}'.split()
    process_report = subprocess.run(cmd, capture_output=True)
    gflops = float(process_report.stdout.split()[3])
    return gflops

def gen_data(mode, n, start_ts, end_ts, nsteps, exp_growth=False, nrep=1):
    """Generate performance data over range of tile sizes
    
    mode:              one if 'ijk', 'jik', 'ikj', 'kij', 'jki', or 'kji'
    n:                 matrix dimension - matrices are n x n
    start_ts, end_ts:  initial and final tile sizes 
    nsteps:            number of steps between initial and final tile sizes
    exp_growth:        True for exp. growing tile sizes, False for linear
    nrep:              number of samples to collect; best value is reported    
    """

    # Sampling across the interval [start_ts, end_ts] can be done either
    # in linear (equally spaced) steps or in exponentially growing steps.
    if exp_growth:
        factor = (end_ts / start_ts)**(1 / nsteps)
    else:
        step = (end_ts - start_ts) / nsteps

    # Initalize lists to hold tile sizes and corresponding GFlops
    ts = []
    gflops = []

    for i in range(nsteps+1):
        if exp_growth:
            ts.append(int(np.round(start_ts * factor**i)))
        else:
            ts.append(int(np.round(start_ts + i*step)))
        gflop_max = 0.0
        for k in range(nrep):
            # keep largest GFlop/s rate
            gflop_max = max(gflop_max, get_gflops(ts[-1], mode, n))
        gflops.append(gflop_max)
        print(f'{ts[-1]} {gflops[-1]}')
    return ts, gflops

if __name__ == '__main__':
    default_size = 2000
    # set up argument parser
    parser = argparse.ArgumentParser(description='Searches parameter space')
    parser.add_argument('ijk_mode', help='one of ijk, ikj, jik, jki, kij, kji')
    parser.add_argument('start_tile_size', help='initial tile size')
    parser.add_argument('end_tile_size', help='final tile size')
    parser.add_argument('nsteps', help='number of steps')
    parser.add_argument('-n', '--size', dest='n', type=int,
                        default=default_size,
                        help=f'matrices are NxN (default is N={default_size})')
    parser.add_argument('-e', '--exp_growth', action='store_true',
                        help='use exponential steps (default is linear)')
    parser.add_argument('-g', '--graph_data', action='store_true',
                        help='plot GFlop/s vs tile size')
    parser.add_argument('-r', '--repetitions', type=int, default=1,
                        help='Number of trials (fastest trial is reported)')
    args = parser.parse_args()

    # process command line
    mode = args.ijk_mode
    start_ts = int(args.start_tile_size)
    end_ts = int(args.end_tile_size)
    nsteps = int(args.nsteps)
    n = int(args.n)
    exp_growth = args.exp_growth
    graph_data = args.graph_data
    nrep = args.repetitions

    # run trials
    ts, gflops = gen_data(mode, n, start_ts, end_ts, nsteps, exp_growth, nrep)

    ## convert tile sizes from bytes to kilobytes
    #for i in range(len(ts)):
    #    ts[i] = int(np.round(ts[i] / 1024))

    if graph_data:
        # plot results
        plt.plot(ts, gflops, '+-')
        if exp_growth:
            plt.xscale('log')
        #plt.xlabel('Tile Size (kilobytes)')
        plt.xlabel('Tile Size (bytes)')
        plt.ylabel('GFlop/s')
        title = f'{n}x{n} matrix product using {mode} method - '
        title += 'Exponential' if exp_growth else 'Linear' + ' steps'
        plt.title(title)
        plt.show()
