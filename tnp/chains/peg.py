from warnings import warn

import mbuild as mb
from mbuild import clone
from mbuild.lib.moieties import CH2 
from tnp.moieties.one_port.ch3 import CH3
from tnp.moieties.peg import PegMonomer

class PEG(mb.Compound):
    """A polyethylen glycol which may optionally end with a hydrogen or a Port."""
    def __init__(self, n=3, cap_front=True, cap_end=True):
        """Initialize an PEG Compound.

        Args:
            n: Number of carbon + oxygen atoms.
                After deducted for cap front, cap end, and ch2 end
                n will be divided by 3 since there are 3 element in a PEG
            cap_front: Add methyl group to beginning of chain ('down' port).
            cap_end: Add methyl group to end of chain ('up' port).
        """
        if n < 1:
            raise ValueError('n must be 1 or more')
        elif n == 1:
            if cap_front or cap_end:
                warn('Returning CH3 group')
            else:
                warn('Returning CH2 group')
        super(PEG, self).__init__()
        
        if n > 1:
            if cap_front:
                n -= 1
                ch3_front = CH3()
                self.add(ch3_front, 'methyl_front')
            if cap_end:
                n -= 1
                ch3_end = CH3()
                self.add(ch3_end, 'methyl_end')

            # If n is not divisible by 3, add CH2 in place of them
            m = 0
            while n % 3 != 0:
                m += 1
                n -= 1 

            n = n//3
            chain = mb.lib.recipes.Polymer(PegMonomer(), n=n, port_labels=('up', 'down'))
            self.add(chain, 'chain')

            if m != 0:
                ch2_chain = mb.lib.recipes.Polymer(CH2(), n=m, port_labels=('up','down'))
                self.add(ch2_chain, 'ch2_chain')
                mb.force_overlap(self['ch2_chain'], self['ch2_chain']['up'],
                                self['chain']['down'])
                if cap_end:
                    mb.force_overlap(self['methyl_end'], self['methyl_end']['down'],
                                    self['ch2_chain']['down'])
                else:
                    self.add(ch2_chain['down'], 'down', containment=False)
            else: 
                if cap_end:
                    mb.force_overlap(self['methyl_end'], self['methyl_end']['down'],
                                    self['chain']['down'])
                else:
                    self.add(chain['down'], 'down', containment=False)
            if cap_front:
                mb.force_overlap(self['methyl_front'], self['methyl_front']['down'],
                                self['chain']['up'])
            else:
                self.add(chain['up'], 'up', containment=False)

        else:
            if cap_end or cap_front:
                ch3 = CH3()
                self.add(ch3, 'methyl')
                self.add(ch3['down'], 'down', containment=False)
            else:
                peg = PegMonomer()
                self.add(peg, 'PEG')
                self.add(peg['up'], 'up', containment=False)
                self.add(peg['down'], 'down', containment=False)

