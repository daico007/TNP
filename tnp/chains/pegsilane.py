import mbuild as mb
from mbuild.lib.moieties import Silane

from tnp.chains.peg import PEG
from tnp.moieties.one_port.ch3 import CH3
class PEGSilane(mb.Compound):
    """A terminal-functionalized pegsilane chain.

    An pegsilane chain featuring a user-specified functional group at one
    terminus and a silane group (featuring an open port to attach to a surface)
    at the other terminus.

    Parameters
    ----------
    n : int
        Length of the chain (number of carbons)
    ch2_cap : True
        Optional CH2 cap for the chain
    """
    def __init__(self, n, ch2_cap=False):
        super(PEGSilane, self).__init__()

        tgroup = CH3() 

        peg = PEG(n, cap_front=False, cap_end=False)
        self.add(peg, 'peg')
        self.add(tgroup, 'terminal_group')
        mb.force_overlap(self['peg'], self['peg']['up'], 
                         self['terminal_group']['down'])
        silane = Silane()
        self.add(silane, 'silane')
        mb.force_overlap(self['silane'], self['silane']['up'], self['peg']['down'])

        self.add(silane['down'], 'down', containment=False)

if __name__ == '__main__':
    chain = PEGSilane(n=8, terminal_group = 'methyl')
    chain.save('pegsilane.pdb', overwrite=True)

