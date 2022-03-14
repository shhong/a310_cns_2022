
import numpy as np
from neuron import h, gui


from cell_template import Cell
class ConnorStevens(Cell):

    def __init__(self, i):
        # We initialize the cell as in the `Cell` class.
        super().__init__()
        
        # Then, we can add additional operations,
        # such as giving each cell its own "kickstart" stimulator,
        self.stim = []
        self.stim.append(h.NetStimFD(self.soma(0.5)))
        self.stim[-1].interval = 3
        self.stim[-1].noise = 1
        self.stim[-1].start = 0
        self.stim[-1].duration = 7
        self.stim[-1].seed(i+1223)

        self.stim.append(h.NetCon(self.stim[0], self.synlist[0]))
        self.stim[-1].weight[0] = 5e-3

    def create_sections(self):
        """create a soma"""
        self.soma = h.Section(name="soma", cell=self)

    def build_topology(self):
        pass # single compartment

    def build_subsets(self):
        pass # single compartment

    def define_geometry(self):
        self.soma.L = 15
        self.soma.diam = 15

    def define_biophysics(self):
        h.v_init = -65
        self.soma.insert('conste')

    def create_synapses(self):
        self.synlist.append(h.Exp2Syn(self.soma(0.5))) # Excitatory
        self.synlist[-1].tau1 = 0.1
        self.synlist[-1].tau2 = 2.0
        self.synlist[-1].e = 0

        self.synlist.append(h.Exp2Syn(self.soma(0.5))) # Inhibitory
        self.synlist[-1].e = -75
        self.synlist[-1].tau1 = 0.1
        self.synlist[-1].tau2 = 2.0


Ncells = 1000                    # Number of cells
Nexc = int((Ncells/5)*4)         # Excitatory cells = 80%
Ninh = int(Ncells/5)             # Inhibitory cells = 20%
nexcpre = int(Nexc*0.1)
ninhpre = int(Ninh*0.1)


from net_manager import ParallelNetManager # Note we are using the parallel version
pnm = ParallelNetManager(Ncells)

for i in range(Ncells):
    pnm.register_cell(i, ConnorStevens(i))

    
gexc = 3e-3      # g_exc = 3 nS
ginh = 5.*gexc   # g_inh/g_exc = 5


for i in pnm.gidlist:
    exc_pre = np.random.randint(0, Nexc, nexcpre)     # Randomly choose 10% of the exc cells
    inh_pre = np.random.randint(Nexc, Ncells, ninhpre) # Randomly choose 10% of the inh cells

    for k in exc_pre:
        if i!=k:  # No self-connection
            pnm.nc_append(k, i, 0, gexc, 0.1, thresh=-10)

    for k in inh_pre:
        if i!=k:  # No self-connection
            pnm.nc_append(k, i, 1, ginh, 0.1, thresh=-10)

pnm.want_all_spikes(thresh=-10)
pnm.set_output_filename('spikes_1.dat')

h.tstop = 250
h.init()
pnm.run()

pnm.message("Simulation complete!")

h.quit()
