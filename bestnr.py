import numpy as np


def bestnr(coherent_snr, sngl_snr, bank_chisq, bank_dof, auto_chisq, auto_dof,
           power_chisq, power_dof, snr_threshold=6, sngl_thresh=4.0, chisq_threshold=None):
    """
    An attempt to replicate the BestNR algorithm from coh_ptf into pycbc_multi_inspiral.
    Rather than use coherent data, this is a test using the single detector chisqs and snrs.
    Some stuff was ignored, I'm just looking to see if the results of the bank and auto chisqs are good.
    Use order H1, L1

    Parameters
    ----------
    coherent_snr: just the number
    sngl_snr: array per ifo
    bank_chisq: array per ifo
    bank_dof: integer
    auto_chisq: array per ifo
    auto_dof: integer
    power_chisq: array per ifo
    power_dof: integer
    snr_threshold: number
    sngl_thresh: number
    chisq_threshold: number

    Returns
    -------
    bestnr: tuple of bestnr per ifo
    """

    # Some variable definitions
    index = 4.0
    nhigh = 3.0

    # Coherent SNR cut
    if coherent_snr < snr_threshold:
        return 0,0

    if chisq_threshold is None:
        chisq_threshold = snr_threshold

    # Reduced chisqs for new snr calculations
    reduced_bank = bank_chisq / bank_dof
    reduced_auto = auto_chisq / auto_dof
    reduced_power = power_chisq / power_dof

    # New bank SNR calculation
    new_bank = np.ones(len(bank_chisq))
    for i in range(len(new_bank)):
        if reduced_bank[i] > 1:
            new_bank[i] = sngl_snr[i] / ((1+reduced_bank[i]**(index/nhigh))/2)**(1.0/index)
        else:
            new_bank[i] = sngl_snr[i]

    # New auto SNR calculation
    new_auto = np.ones(len(auto_chisq))
    for i in range(len(new_auto)):
        if reduced_auto > 1:
            new_auto[i] = sngl_snr[i] / ((1+reduced_auto[i]**(index/nhigh))/2)**(1.0/index)
        else:
            new_auto[i] = sngl_snr[i]

    # New power SNR calculation (used later)
    new_power = np.ones(len(power_chisq))
    for i in range(len(new_power)):
        if reduced_power[i] > 1:
            new_power[i] = sngl_snr[i] / ((1+reduced_power[i]**(index/nhigh))/2)**(1.0/index)
        else:
            new_power[i] = sngl_snr[i]

    # Now cut based on bank and auto
    for i in range(len(new_bank)):
        if (new_bank[i] < chisq_threshold) or (new_auto[i] < chisq_threshold):
            return 0,0

    # Single det snr cut
    # NOTE: pass single det snr as a tuple
    for snr in sngl_snr:
        if snr < sngl_thresh:
            return 0,0

    # NOTE: I'm ignoring the null SNR cuts for now
    # Which means use a 1 or 2 ifo run

    # And now for the final BestNR, just reweighted by power chisq
    return new_power


if __name__ == '__main__':
    import csv
    import h5py
    file_name = 'BestNR_output.csv'
    in_file = input('Enter the data filename:\n')
    hdf_file = h5py.File(in_file)
    net = hdf_file['network']
    h1 = hdf_file['H1']
    l1 = hdf_file['L1']

    with open(file_name, 'wt') as file:
        writer = csv.writer(file)
        header = ['End time', 'Coherent SNR', 'H1 SNR', 'L1 SNR', 'Power chisq H1', 'Power chisq L1',
                  'Bank chisq H1', 'Bank Chisq L1', 'Auto chisq H1', 'Auto chisq L1', 'BestNR H1', 'BestNR L1']
        end_time = []
        coh_snr = []
        snr_h1 = []
        snr_l1 = []
        power_h1 = []
        power_l1 = []
        bank_h1 = []
        bank_l1 = []
        auto_h1 = []
        auto_l1 = []
        bestnr_h1 = []
        bestnr_l1 = []

        # End time
        for i in net['end_time_gc']:
            end_time.append(i)

        # Coherent SNR
        for i in net['coherent_snr']:
            coh_snr.append(i)

        # SNR H1
        for i in h1['snr']:
            snr_h1.append(i)

        # SNR L1
        for i in l1['snr']:
            snr_l1.append(i)

        # Power chisq
        for i in h1['chisq']:
            power_h1.append(i)

        for i in l1['chisq']:
            power_l1.append(i)

        # Bank chisq
        for i in h1['bank_chisq']:
            bank_h1.append(i)

        for i in l1['bank_chisq']:
            bank_l1.append(i)

        # Auto chisq
        for i in h1['cont_chisq']:
            power_h1.append(i)

        for i in l1['cont_chisq']:
            power_l1.append(i)

        # Now get bestNR
        # get chisq dofs
        power_dof = h1['chisq_dof'][0]
        bank_dof = h1['bank_chisq_dof'][0]
        auto_dof = h1['cont_chisq_dof'][0]
        for i in range(len(end_time)):
            # call the function, get the tuple, and put the thing in the array
            bestnr_tuple = bestnr(coh_snr[i], (snr_h1[i], snr_l1[i]), (bank_h1[i], bank_l1[i]), bank_dof,
                                  (auto_h1[i], auto_l1[i]), auto_dof, (power_h1[i], power_l1[i]), power_dof)

            bestnr_h1.append(bestnr_tuple[0])
            bestnr_l1.append(bestnr_tuple[1])

        # Output into csv file
        writer.writerow(header)
        for i in range(len(end_time)):
            writer.writerow([end_time[i], coh_snr[i], snr_h1[i], snr_l1[i], power_h1[i], power_l1[i],
                             bank_h1[i], bank_l1[i], auto_h1[i], auto_l1[i], bestnr_h1[i], bestnr_l1[i]])
