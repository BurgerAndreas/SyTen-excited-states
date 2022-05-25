import numpy as np
import matplotlib.pyplot as plt
import csv
import os as os


# plot energy in file ./log_folder/fname_energies.csv
def plot_energies(fname='ising32', spin='Half', clean=None):
    #
    log_folder = './data_spin' + spin + '/'
    plot_folder = './plots_spin' + spin + '/'
    #log_folder = './data/'
    #plot_folder = './plots/'

    # readout energies
    log_file = log_folder + fname + '_energy.csv'
    energies_file = np.loadtxt(fname=log_file, delimiter=',', dtype=np.cdouble)
    j = np.real(energies_file[:, 0])
    energies = np.real(energies_file[:, 1:])

    # sort energies, in case GS converged to an excited state
    for cnt, row in enumerate(energies):
        energies[cnt] = np.sort(energies[cnt])

    # merge 'degenerate' states
    if clean is not None:
        energies_reduced = np.zeros([np.shape(energies)[0], int(np.shape(energies)[1]/2)])
        if type(clean) is list:
            for e_num, e in enumerate(clean):
                energies_reduced[:, e_num] = energies[:, e]
        else:
            for e in range(np.shape(energies)[1]):
                if e % 2 == 1:
                    for t in range(np.shape(energies)[0]):
                        # take lower one
                        if clean == 'smaller':
                            if energies[t, e] < energies[t, e-1]:
                                energies_reduced[t, int(e/2)] = energies[t, e]
                            else:
                                energies_reduced[t, int(e/2)] = energies[t, e-1]

                        # take more consistent one
                        elif clean == 'smaller_tangent':
                            tangent1 = energies[t-1, e - 1] - energies[t, e - 1]
                            tangent2 = energies[t-1, e] - energies[t, e]
                            if np.abs(tangent2) < np.abs(tangent1):
                                energies_reduced[t, int(e / 2)] = energies[t, e]
                            else:
                                energies_reduced[t, int(e / 2)] = energies[t, e - 1]

                        elif clean == 'closer_to_tangent':
                            if t == np.shape(energies)[0] - 1:
                                # smaller tangent
                                tangent1 = energies[t - 1, e - 1] - energies[t, e - 1]
                                tangent2 = energies[t - 1, e] - energies[t, e]
                                if np.abs(tangent2) < np.abs(tangent1):
                                    energies_reduced[t, int(e / 2)] = energies[t, e]
                                else:
                                    energies_reduced[t, int(e / 2)] = energies[t, e - 1]
                            else:
                                tangent1 = energies[t+1, e-1] - energies[t-1, e-1]
                                tangent2 = energies[t+1, e] - energies[t-1, e]
                                # distance from whats predicted by tangent
                                d1 = np.abs((energies[t-1, e-1] + tangent1/2) - energies[t, e-1])
                                d2 = np.abs((energies[t-1, e] + tangent1/2) - energies[t, e])
                                if d2 < d1:
                                    energies_reduced[t, int(e / 2)] = energies[t, e]
                                else:
                                    energies_reduced[t, int(e / 2)] = energies[t, e - 1]

        energies = energies_reduced

    # plot energies
    for e in range(np.shape(energies)[1]):
        plt.plot(j, energies[:, e], marker='x', lw=1, ms=4, label='E' + str(e))
    # plt.grid(axis='both')
    plt.title('Energies | ' + fname + ' spin' + spin)
    plt.ylabel('energy')
    plt.xlabel('J')
    plt.legend(loc="upper right")
    plt.savefig(plot_folder + fname + '_energies.png', dpi=400)
    plt.show()
    plt.close()
    plt.cla()
    plt.clf()

    # plot energy differences
    for e in range(np.shape(energies)[1] - 1):
        diff = np.abs(energies[:, 0] - energies[:, e + 1])
        plt.plot(j, diff, marker='x', lw=1, ms=4, label='|E0-E' + str(e+1) + '|')
    plt.grid(axis='both')
    plt.title('Energy gaps | ' + fname + ' spin' + spin)
    plt.ylabel('energy gap')
    plt.xlabel('J')
    plt.legend(loc="upper right")
    plt.savefig(plot_folder + fname + '_gap.png', dpi=400)
    plt.show()
    plt.close()
    plt.cla()
    plt.clf()

    # read timings
    log_file = log_folder + '/' + fname + '_time.csv'
    time_file = np.loadtxt(fname=log_file, delimiter=',', dtype=np.cdouble)
    j_time = np.real(time_file[:, 0])
    time = np.real(time_file[:, 1:])

    # plot timings
    for e in range(np.shape(time)[1]):
        plt.plot(j_time, time[:, e], marker='x', lw=1, ms=4, label='E' + str(e))
    # plt.grid(axis='both')
    plt.title('Timings | ' + fname + ' spin' + spin)
    plt.ylabel('time [s]')
    plt.xlabel('J')
    plt.legend(loc="upper right")
    plt.savefig(plot_folder + fname + '_time.png', dpi=400)
    #plt.show()
    plt.close()
    plt.cla()
    plt.clf()

    # save in cleaned format
    j = np.expand_dims(j, axis=0)
    j_energies = np.append(j.T, energies, axis=1)
    j_time = np.expand_dims(j_time, axis=0)
    j_timing = np.append(j_time.T, time, axis=1)

    clean_folder = './data_spin' + spin + '/sorted/' + fname
    np.savetxt(fname=clean_folder + '_energy.csv', X=j_energies, delimiter=',')
    np.savetxt(fname=clean_folder + '_time.csv', X=j_timing, delimiter=',')
    # np.loadtxt(fname='./ising32_energy.csv', delimiter=',', dtype=np.double)


# ising8, heisenberg16
def plot_all():
    for model in ['ising']:
        for sites in [4, 8, 16, 32]:
            for spin in ['Half']:
                plot_energies(fname=model+str(sites), spin=spin)

plot_energies(fname='ising4', spin='Half', clean=[0, 3, 5])
#plot_energies(fname='ising8', spin='Half', clean=[0, 3, 5])
#plot_energies(fname='ising16', spin='Half', clean=[0, 3, 4])
#plot_energies(fname='ising32', spin='Half', clean='smaller')
#plot_energies(fname='ising64', spin='Half', clean='smaller')
#plot_energies(fname='heisenberg32')




