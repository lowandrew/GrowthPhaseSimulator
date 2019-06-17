# Growth Phase Simulator

We tend to sequence our stuff in growth phase, and so end up having lots more sequence close 
to the origin of replication than far away - I'm not aware of any simulators that take this into
account, so this acts as a wrapper around ART (tested with version 2016-06-05 MountRainier) that
will automatically find the origin of replication in _E. coli_ and then simulate reads such that there's
more sequencing coverage around the origin than elsewhere.

## Installation

Clone this repository and use the `growth_phase_simulator.py` script. You'll need to have `art_illumina` 
accessible on your $PATH, as well as `blastn` and `makeblastdb`. In terms of python packages, all that you 
need is biopython.

## Usage

Coming soon... 
