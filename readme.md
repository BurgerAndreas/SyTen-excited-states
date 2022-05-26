# SyTen excited states

Finding excited states using the SyTen library <br />
Examples: 1d Ising and 1d Heisenberg model


SyTen library <br />
by Claudius Hubig, Felix Lachenmaier, Nils-Oliver Linden, Teresa Reinhard, Leo Stenzel, Andreas Swoboda, Martin Grundner and Sam Mardazad <br />
https://syten.eu/

Code in large parts by Sebastian Packel and Sam Mardazad


### Datasets
Ising spinHalf <br />
hamilton.append( - 2 * J_I * lattice.get("sz",i)*lattice.get("sz",i+1) )  # zz <br />
hamilton.append( - h_I  *  lattice.get("sx",i) )  # x

Ising spinHalf2 <br />
hamilton.append( - 4 * J_I * lattice.get("sz",i)*lattice.get("sz",i+1) )  # zz <br />
hamilton.append( - 2 * h_I  *  lattice.get("sx",i) )  # x

Ising spinHalf3 <br />
hamilton.append( - 4 * J_I * lattice.get("sz",i)*lattice.get("sz",i+1) )  # zz <br />
hamilton.append( - 2 * np.sqrt(2) * h_I  *  lattice.get("sx",i) )  # x

### Timings
Intel(R) Xeon(R) CPU E5-2630 v4 at 2.20GHz <br />
