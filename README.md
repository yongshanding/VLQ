# Fault-tolerant non-clifford gates using error detection

Credit: Source code repository is largely built based on the software from
[Virtualized Logical Qubits: A 2.5D Architecture for Error-Corrected Quantum Computing](https://arxiv.org/abs/2009.01982)
published in the proceedings of MICRO '20, the 53rd IEEE/ACM International Symposium on Microarchitecture, October 2020. Forked from [cduck/VLQ](https://github.com/cduck/VLQ)



## Install

1. Clone this repository
    ```bash
    git clone https://github.com/yongshanding/vlq
    cd vlq
    ```

2. Install Julia (tested with 1.4.2): [julialang.org/downloads](https://julialang.org/downloads/)

3. Set path environment. In `~/.bashrc`:
    ```bash
    export PATH="$PATH:/path/to/<Julia directory>/bin"
    ```
    
3. Install required Julia packages (run from the `vlq/` directory) from the julia REPL:
    ```julia
    Pkg.update()
    ENV["PYTHON"]=""; Pkg.build("PyCall")
    ] activate .; instantiate
    ```

4. (Optional) Install packages globally from julia REPL:
    ```julia
    add LightGraphs BlossomV ChpSim OrderedCollections PyPlot
    ```

## Usage

Run the following Julia code in a REPL (start one with `julia --project=.`)

```julia
# Setup
include("src/make_plots.jl"); using .MakePlots
MakePlots.setup(num_workers=16)  # Setup 16 worker processes

# Run the simulations
job_id = "rz_plots0.001"
samples = 2000*100 # How many samples, ~1.15 CPU-hours per 100000 samples
dists = [3:2:9...] # code distances
plots = [1:3...] # types of plots
MakePlots.dist_calc_all(job_id, dists, samples, plots) # Start computing with workers in background
MakePlots.fetch_finished(job_id, 2) # Obtain results as they are computed 
```

