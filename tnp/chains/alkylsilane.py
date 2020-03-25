import mbuild as mb
from mbuild.lib.recipes.alkane import Alkane
from mbuild.lib.moieties.silane import Silane
from tnp.moieties.one_port.ch3 import CH3


class Alkylsilane(mb.Compound):
    """A terminal-functionalized alkylsilane chain.

    An alkylsilane chain featuring a user-specified functional group at one
    terminus and a silane group (featuring an open port to attach to a surface)
    at the other terminus.

    Parameters
    ----------
    n : int
        Length of the chain (number of carbons)
    """
    def __init__(self, n):
        super(Alkylsilane, self).__init__()

        tgroup = CH3()

        alkane = Alkane(n, cap_front=False, cap_end=False)
        self.add(alkane, 'alkane')
        self.add(tgroup, 'terminal_group')
        mb.force_overlap(self['alkane'], self['alkane']['up'], 
                         self['terminal_group']['down'])
        silane = Silane()
        self.add(silane, 'silane')
        mb.force_overlap(self['silane'], self['silane']['up'], self['alkane']['down'])

        self.add(silane['down'], 'down', containment=False)

if __name__ == "__main__":
    chain = Alkylsilane(n=8, terminal_group='cyclopropyl')
    chain.save('alkylsilane.mol2', overwrite=True)

