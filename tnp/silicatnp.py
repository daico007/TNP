from __future__ import division

import mbuild as mb
import numpy as np

from tnp.chains import Alkane, Alkylsilane, PEG, PEGSilane
from tnp.nanoparticles import SilicaSphere
from tnp.nanoparticles import SilicaNP

class SilicaTNP(mb.Compound):
    """A tethered (spherical) nanoparticle

    A tethered nanoparticle, with the choice of
    Alkylsilane and PEGSilane chains.

    Parameters
    ----------
    radius : float, optional, default=4
        Radius of the nanoparticles
    chain_prototype : str, optional, default='peg'
        Can choose from 'alkane', 'alkylsilane', 'peg',
        and 'pegsilane'
    chain_density : float, optional, default=0.5
        Density of the tethered chain (from 0 to 1)
    chain_length : int, optional, default=17
        Chain_length of the tethered chains
    """
    def __init__(self, radius=4, chain_prototype='pegsilane',
            chain_density=0.5, chain_length=17):
        super(SilicaTNP, self).__init__()

        surface_area = 4.0 * np.pi * radius ** 2.0
        n_chains = int(chain_density * surface_area)
        pattern = mb.SpherePattern(n_chains)
        silica_radius = 0.201615

        core_particles = ['Si', 'O']
        shell = SilicaSphere(n=n_chains, radius=radius)
        shell.translate_to((0,0,0))
        core = SilicaNP(radius=radius)
        core.translate_to((0,0,0))
        self.add(shell, 'shell')
        self.add(core, 'core')

        if chain_prototype == 'alkane':
            chain_prototype = Alkane(n=chain_length,
                                 cap_end=False)
        elif chain_prototype == 'alkylsilane':
            chain_prototype = Alkylsilane(n=chain_length)
        elif chain_prototype == 'peg':
            chain_prototype = PEG(n=chain_length,
                                cap_end=False)
        elif chain_prototype == 'pegsilane':
            chain_prototype = PEGSilane(n=chain_length)
        else:
            raise MbuildError('chain prototype not supported')

        # This radius is for an alkane chain,
        # PEG probably close enough to this
        chain_bead_radius = 0.4582 / 2
        port_separation = np.linalg.norm(chain_prototype['down'].pos - \
                                chain_prototype['down'].anchor.pos)
        pattern.scale(radius + silica_radius + chain_bead_radius - \
                      port_separation)


        chains, _ = pattern.apply_to_compound(guest=chain_prototype,
            guest_port_name='down', host=self['shell'])
        self.add(chains)

        """
        self.label_rigid_bodies(rigid_particles=core_particles)

        for bond in self.bonds():
            if(bond[0].name == 'SilicaNP' or
                    bond[1].name == 'SilicaNP'):
                if bond[0].name == 'SilicaNP':
                    bond[1].rigid_id = 0
                if bond[1].name == 'SilicaNP':
                    bond[0].rigid_id = 0
                self.remove_bond(bond)
        """

